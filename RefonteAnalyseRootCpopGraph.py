"""
Script allowing to convert .root raw data of Geant4 in data of interest

Usage
======
    Run script to open graphical interface

Returns
=======
    doses to cell nucleus and cytoplams
    cell survivals
    cross-fire information
"""

# from geometry_informations import *
import geometry_informations
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas
from pyexcelerate import Workbook
import scipy.integrate
import scipy.interpolate as interpolate
import sys
import time
import tkinter
import tkinter.ttk
import uproot
import warnings
from xml.dom import minidom

warnings.filterwarnings("error")

KEV_IN_J = 1.60218 * 1e-16
WATER_DENSITY = 1e-3  # 10-3 kg/cm³
UNIT_COEFFICIENT = (KEV_IN_J / (WATER_DENSITY * (1e-4)))

print(UNIT_COEFFICIENT)

SIG0 = [49*np.pi, 24.01*np.pi, 34.81*np.pi] #From Mario calculations
A_CST = 0.1602 #Gy μm3 keV−1
energies_valid_for_alpha_beta_approximation = np.arange(200,90001)

################ Fit parameters used in article TCP RIV-alpha ###################
# Y0 = [0.06072969, 0.02562553, 0.03934994] #HSG, V79, CHO-K1
# A = [-0.18385472, -0.10426184, -0.11163773]
# W = [3.05093045, 2.87758559, 3.20398251]
# XC = [0.46545609, 0.38084839, 0.48452192]

################ Fit parameters from 2022/12/16 ################### TO DO : evaluate differences with previous parameters
Y0 = [0.06486164, 0.02722410, 0.04221387] #HSG, V79, CHO-K1
A = [-0.26336407, -0.11801719, -0.19357751]
W = [3.39940424, 2.97713123, 3.90866411]
XC = [-0.00863166, 0.23348883, -0.25238105]

BETAG = [0.0961, 0.0405, 0.0625]  # constante de Monini et al. 2019


Emax=8000 #Energie max des ions Hélium émis, en keV

radius_cell_line = 7 * 1e-4  # Rayon du noyau de la lignée HSG
surface_centerslice_cell_line = math.pi * radius_cell_line ** 2
length_of_cylinderslice_cell = 1

bins = 200
START_TIME = time.time()

np.set_printoptions(threshold=sys.maxsize)
#np.set_printoptions(threshold = False)

global labeling_percentage, Labeling_Percentage_Entry, cell_compartment_radionuclide_Txt, distrib_name_cb

def column(matrix, i):
    return [row[i] for row in matrix]

def conversion_energy_in_let(data_base, energy):
    """
    Returns a function that converts an input energy into the corresponding LET from a given data base

    Input
    -------
    data_base in string format
    energy in keV

    """

    TABLES_CONVERSION_ENERGY_IN_LET = pandas.read_excel("E_TEL/conversion_tables_" + data_base +".xlsx").to_records()
    print("E_TEL/conversion_tables_" + data_base +".xlsx")
    ENERGY_LIST = TABLES_CONVERSION_ENERGY_IN_LET['E(keV)']
    CORRESPONDING_LET_LIST = TABLES_CONVERSION_ENERGY_IN_LET['LET(keV/um)']
    continuous_function_to_convert_energy_in_let = interpolate.interp1d(ENERGY_LIST, CORRESPONDING_LET_LIST, fill_value="extrapolate", kind= "linear")

    return continuous_function_to_convert_energy_in_let(energy)

def beta_nanox(E,type_cell): #E in MeV/nucleon
    return Y0[type_cell] + (A[type_cell]/(2*np.pi))*(W[type_cell]/(((E-XC[type_cell])**2)+(W[type_cell]**2)/4))

def alpha_nanox(E,type_cell):
    conv_LET_E_SRIM = conversion_energy_in_let("SRIM", E)
    b = BETAG[type_cell] * (((A_CST * conv_LET_E_SRIM * 0.8)/SIG0[type_cell]) ** 2)
    return ((SIG0[type_cell] +((A_CST * conv_LET_E_SRIM * (b-1))*np.sqrt(beta_nanox(E/4000,type_cell)/(b+(b*b)/2))))/(A_CST * conv_LET_E_SRIM))

def dn1_dE_continous():
    conversion_energy_in_let_srim_alpha_beta_approximation_range = \
        conversion_energy_in_let("SRIM", energies_valid_for_alpha_beta_approximation)
    conversion_energy_in_let_g4_alpha_beta_approximation_range = \
        conversion_energy_in_let("G4", energies_valid_for_alpha_beta_approximation)
    dn1_dE = -np.log(1 - alpha_nanox(energies_valid_for_alpha_beta_approximation,0) \
                     * UNIT_COEFFICIENT*conversion_energy_in_let_srim_alpha_beta_approximation_range \
                     / surface_centerslice_cell_line) \
                     / (length_of_cylinderslice_cell * conversion_energy_in_let_g4_alpha_beta_approximation_range) #calculation of number of lethal events per keV, via Mario's approximations
    dn1_dE_interpolated = interpolate.interp1d(energies_valid_for_alpha_beta_approximation, dn1_dE, fill_value="extrapolate", kind="linear")
    return dn1_dE_interpolated

def number_of_lethal_events_for_alpha_traversals(dn1_dE_function):
    """
    Returns the function that converts an energy E into the cumulated number of lethal damage from 0 to E
    """
    energie_table_binned = np.linspace(0, Emax, num=bins)
    f_He_cumulative_int = scipy.integrate.cumtrapz(dn1_dE_function(energie_table_binned),
                                             energie_table_binned, initial=0)
    n1 = interpolate.interp1d(energie_table_binned, f_He_cumulative_int, fill_value="extrapolate",
                              kind="linear")  # fonction primitive continue en fonction de E
    return n1

def determine_cells_in_2_spheroid_zones(positions_x, positions_y, positions_z, radius_zone_1, radius_zone_2, nb_cells):
    """
    Returns the number of cells in zone 1 and 2 for chosen radii (in µm),
    and an array, sorted by cell id, with the zones where the cells are
    """
    positions_cell = np.sqrt(positions_x ** 2 + positions_y ** 2 + positions_z ** 2)
    zone_cell = np.zeros(nb_cells)
    nb_cell_zone_1 = 0
    nb_cell_zone_2 = 0
    index_list = np.arange(0,nb_cells)
    for index in index_list:
        if (positions_cell[index] < radius_zone_1):
            zone_cell[index] = 1
            nb_cell_zone_1 += 1
        elif (positions_cell[index] < radius_zone_2):
            zone_cell[index] = 2
            nb_cell_zone_2 += 1
    return(zone_cell, nb_cell_zone_1, nb_cell_zone_2)

def todo_histo():
    histo_nb_noy_par_p = 0 # mettre 1 affiche l'histogramme du nombre de noyaux traversés par les particules
    histo_edep_noy_par_p = 0 # mettre 1 affiche l'histogramme de l'énergie moyenne déposée par particule dans un noyau quand elle y rentre

    if (histo_nb_noy_par_p==1) :
        _, _, patches = plt.hist(nb_nucl_traversées_par_la_particule_tab_sur_toutes_simus_np, bins=100, edgecolor='white')

        plt.xticks(fontsize=13.5)
        plt.yticks(fontsize=13.5)
        plt.xlabel('Nb of nuclei crossed by a particle',fontsize=15,fontname="Liberation Sans",fontweight='bold')
        plt.ylabel('Occurence',fontsize=15,fontname="Liberation Sans",fontweight='bold')
        #plt.title('Histogram of : Nb of nuclei crossed by a particle')
        plt.grid(True)
        plt.show()

    if (histo_edep_noy_par_p==1) :
        Edep_dans_noy_par_particule_np_histo=np.resize(Edep_dans_noy_par_particule_np,(1,len(Edep_dans_noy_par_particule_np)*len(Edep_dans_noy_par_particule_np[0])))
        _, _, patches = plt.hist(Edep_dans_noy_par_particule_np_histo[0], bins=1000, edgecolor='black')

        plt.xticks(fontsize=13.5)
        plt.yticks(fontsize=13.5)
        plt.xlabel('Edep by a particle in cell nucleus (keV)',fontsize=15,fontname="Liberation Sans",fontweight='bold')
        plt.ylabel('Occurence',fontsize=15,fontname="Liberation Sans",fontweight='bold')
        # plt.title('Histogram of : Edep by a particle in cell nucleus')
        plt.grid(True)
        plt.show()

def if_internalization_study():
    if labeling_percentage.winfo_exists():
        labeling_percentage.destroy()
    if Labeling_Percentage_Entry.winfo_exists():
        Labeling_Percentage_Entry.destroy()
    cell_compartment_radionuclide_Txt = tkinter.Label(window, text="Intra cellular distribution name :", fg='blue')
    cell_compartment_radionuclide_Txt.place(x=100, y=100)
    selected_distrib_name = tkinter.StringVar()
    distrib_name_cb = tkinter.ttk.Combobox(window, width=35, textvariable=selected_distrib_name)
    distrib_name_cb['values'] = ['Membrane', 'Cytoplasm', 'Homogeneous', 'Nucleus']
    distrib_name_cb.place(x=400, y=100)

def if_labeling_study() :
    if cell_compartment_radionuclide_Txt.winfo_exists():
        cell_compartment_radionuclide_Txt.destroy()
    if distrib_name_cb.winfo_exists():
        distrib_name_cb.destroy()
    labeling_percentage = tkinter.Label(window, text="Labeling percentage : ", fg='blue')
    labeling_percentage.place(x=100, y=100)
    Labeling_Percentage_Entry = tkinter.Entry(window, width=35)
    Labeling_Percentage_Entry.place(x=400, y=100)

def id_deletion_of_root_outputs_with_errors():
    """
    Returns
    ======
        ids of simulations that returned a root output without errors when they were opened, in list format
        ids of simulations that returned a root output without errors when they were opened, in numpy format
        the number of errors encountered
    """
    nb_opened_files_for_one_simulation = indexe_of_root_output = 0
    index_complete_simulation = 1
    hfin = nb_files_for_one_complete_simulation
    indexes_root_files_without_errors_temp, indexes_root_files_without_errors = ([] for nb_arrays in range(2))

    while (indexe_of_root_output < (nb_complete_simulations * nb_files_for_one_complete_simulation)):
        try:
            if ((((indexe_of_root_output - 1) % nb_files_for_one_complete_simulation) == 0) & (indexe_of_root_output != 1)):
                index_complete_simulation += 1

            root_file_name = "Root/" + "outputMultiCellulaire/" + dossier_root + nom_fichier_root +\
                            f"{indexe_of_root_output}" + "_t0" + ".root"
            f1 = uproot.open(root_file_name)
            d1 = f1['cell']['nameParticle'].array(library="np")
            indexes_root_files_without_errors_temp.append(indexe_of_root_output)
            nb_opened_files_for_one_simulation += 1
            indexe_of_root_output += 1

            if (nb_opened_files_for_one_simulation == nb_files_for_one_complete_simulation):
                indexes_root_files_without_errors.append(indexes_root_files_without_errors_temp)
                indexes_root_files_without_errors_temp = []
                nb_opened_files_for_one_simulation = 0

        except:
            indexe_of_root_output = (index_complete_simulation * nb_files_for_one_complete_simulation) + 1
            index_complete_simulation += 1
            nb_opened_files_for_one_simulation = 0
            indexes_root_files_without_errors_temp = []

    indexes_root_files_without_errors_np = np.asarray(indexes_root_files_without_errors)
    indexes_root_files_without_errors_np = np.resize(indexes_root_files_without_errors_np,
                               (1,(len(indexes_root_files_without_errors_np[0]) * len(indexes_root_files_without_errors_np))))
    indexes_root_files_without_errors_np = np.sort(indexes_root_files_without_errors_np[0])
    nb_files_with_errors = nb_complete_simulations - len(indexes_root_files_without_errors)

    return(indexes_root_files_without_errors, indexes_root_files_without_errors_np, nb_files_with_errors)

def graphic_window():
    global window, radiovalue_study_type, distrib_name_cb, geom_cb, radionuclide_name_entry, nb_simulations_entry, \
           diffusion_list, number_particles_per_cell_list, cell_line_cb

    window = tkinter.Tk()
    window.geometry("1000x700")

    StudyType_Txt = tkinter.Label(window, text = "Type of study :", fg='red')
    StudyType_Txt.place (x=100, y=50)

    labeling_percentage = tkinter.Label(window, text="Labeling percentage : ", fg='blue')
    Labeling_Percentage_Entry = tkinter.Entry(window, width=35)

    radiovalue_study_type=tkinter.IntVar()
    radiovalue_study_type.set(0)
    r1=tkinter.Radiobutton(window, text="Internalization", variable=radiovalue_study_type,value=0,command=if_internalization_study)
    r2=tkinter.Radiobutton(window, text="Labeling", variable=radiovalue_study_type,value=1, command=if_labeling_study)
    r1.place(x=390,y=50)
    r2.place(x=590, y=50)

    cell_compartment_radionuclide_Txt = tkinter.Label(window, text="Intra cellular distribution name :", fg='blue')
    cell_compartment_radionuclide_Txt.place(x=100, y=100)
    selected_distrib_name = tkinter.StringVar()
    distrib_name_cb = tkinter.ttk.Combobox(window, width=35 , textvariable=selected_distrib_name)
    distrib_name_cb['values'] = ['Membrane', 'Cytoplasm', 'Homogeneous', 'Nucleus']
    distrib_name_cb.current(0)
    distrib_name_cb.place(x=400, y=100)

    geom_name_txt = tkinter.Label(window, text = "Geometry name :", fg='blue')
    geom_name_txt.place (x=100, y=150)
    selected_geom = tkinter.StringVar()
    geom_cb = tkinter.ttk.Combobox(window, width=35 , textvariable=selected_geom)
    geom_cb['values'] = ["30µmRadius Spheroid, 75 % cell packing", "50µmRadius Spheroid, 75 % cell packing", "70µmRadius Spheroid, 75 % cell packing", "160µmRadius Spheroid, 75 % cell packing" ,"95µmRadius Spheroid, 25 % cell packing", "95µmRadius Spheroid, 50 % cell packing", "95µmRadius Spheroid, 75 % cell packing", "95µmRadius Spheroid, 75 % cell packing 2", "100µmRadius Spheroid, 40 % cell packing"]
    geom_cb.place(x=400, y=150)

    radionuclide_name_txt = tkinter.Label(window, text = "Radionuclide used :", fg='blue')
    radionuclide_name_txt.place (x=100, y=200)
    radionuclide_name_entry = tkinter.Entry(window, width=35)
    radionuclide_name_entry.insert(tkinter.END, "At211")
    radionuclide_name_entry.place(x=400, y=200)

    nb_simulations_txt = tkinter.Label(window, text = "Number of simulations to analyse :", fg='blue')
    nb_simulations_txt.place (x=100, y=250)
    nb_simulations_entry = tkinter.Entry(window, width=35)
    nb_simulations_entry.insert(tkinter.END, "20")
    nb_simulations_entry.place(x=400, y=250)

    diffusion_txt = tkinter.Label(window, text = "Daughter diffusion :", fg='blue')
    diffusion_txt.place (x=100, y=300)
    diffusion_choice = tkinter.StringVar()
    diffusion_list = tkinter.ttk.Combobox(window, width=35, textvariable=diffusion_choice)
    diffusion_list['values'] = ["Yes", "No"]
    diffusion_list.current(1)
    diffusion_list.place(x=400, y=300)

    number_particles_per_cell_txt = tkinter.Label(window, text = "Number of alpha particles per cell :", fg='blue')
    number_particles_per_cell_txt.place(x=100, y=350)
    number_particles_per_cell_choice = tkinter.StringVar()
    number_particles_per_cell_list = tkinter.ttk.Combobox(window, width=35, textvariable=number_particles_per_cell_choice)
    number_particles_per_cell_list['values'] = ["1", "2", "3", "4", "5" ,"6", "7", "8", "9" ,"10", "42"]
    number_particles_per_cell_list.current(0)
    number_particles_per_cell_list.place(x=400, y=350)

    cell_line_txt = tkinter.Label(window, text="Cell line :", fg='blue')
    cell_line_txt.place(x=100, y=400)
    selected_cell_line = tkinter.StringVar()
    cell_line_cb = tkinter.ttk.Combobox(window, width=35 , textvariable=selected_cell_line)
    cell_line_cb['values'] = ['HSG', 'V79', 'CHO-K1']
    cell_line_cb.current(0)
    cell_line_cb.place(x=400, y=400)


    Validate = tkinter.Button(window, text = "Validate", command = add_new_buttons_to_graphic_window)
    Validate.place(x=480, y=450)

    window.mainloop()

def create_folder_for_output_analysis_files():
    global dossier_root, index_of_first_root_output, nom_dossier_pour_excel_analyse
    dossier_root = study_type_folder_name + "/" + available_data_name_file[data_cb.current()] + "/"
    index_of_first_root_output = 0 #Works only if the indexes of root files start at 0
    print("dossier_root : ", dossier_root)
    nom_dossier_pour_excel_analyse = available_data_date[data_cb.current()] + "_" + "_" + str(spheroid_compaction) +\
                                     "CP_" + str(r_sph) + "um_" + rn_name + "_diff" +\
                                     bool_diff[diffusion_list.current()] + "_" +\
                                     str(nb_particles_per_cell[number_particles_per_cell_list.current()]) + "ppc" + \
                                     "_" + cell_line_cb.get()
    try:
        os.makedirs(os.path.join("AnalysisResults"))
    except FileExistsError:
        print("Folder AnalysisResults already created")
    try:
        os.makedirs(os.path.join("AnalysisResults/" + study_type_folder_name))
    except:
        print("Folder AnalysisResults/study_type_folder_name already created")
    try:
        os.makedirs(os.path.join("AnalysisResults/" + study_type_folder_name, nom_dossier_pour_excel_analyse))
    except:
        print("Folder for the required analysis results already created")
    return None

def main():
    create_folder_for_output_analysis_files()

    print("selected date : ",available_data_date[data_cb.current()])


    ######################## Conversion des alpha en dn1/dE ########################################

    dn1_dE_continous_pre_calculated = dn1_dE_continous()

    ##################### Gestion des ID de CPOP ##################################################

    txt_id_deleted_cells = "Cpop_Deleted_Cells_ID_Txt/" + "IDCell_" + nom_config + ".txt"

    real_id_cells, test_file_not_empty, deleted_id_txt = geometry_informations.cpop_real_cell_id_determination(txt_id_deleted_cells, nb_cellules_xml)

    nb_cellules_reel = len(real_id_cells)

    print("nb_cellules_reel : ", nb_cellules_reel)

    Perfect_ID_Cells = np.arange(0,nb_cellules_reel)


    ###################### Lecture Geométrie #####################################

    ###### Masses #######

    txt_cells_masses="Cpop_Masse_Txt/" + "MassesCell_" + nom_config + ".txt"

    masses_cytoplasms, masses_nuclei, masses_cells = geometry_informations.masses_cells_reading(txt_cells_masses)

    r_tum = float(r_sph) * 10**(-6) #in meters
    masse_tum=((4/3)*np.pi*r_tum**3)*1000 #in kg

    print("cell_packing = ", sum(masses_cells)/masse_tum)
    print("masse_tum = ", masse_tum)

    ###### Positions ######

    positions_x, positions_y, positions_z = geometry_informations.positions_cells_reading(xml_geom, real_id_cells)

    zone_cell, nb_cell_zone_1, nb_cell_zone_2 = determine_cells_in_2_spheroid_zones(positions_x,
                                                    positions_y, positions_z,
                                                    radius_zone_1 = 50, radius_zone_2 = 95,
                                                    nb_cells = nb_cellules_reel)

    print("nb_cell_zone_1", nb_cell_zone_1)
    print("nb_cell_zone_2", nb_cell_zone_2)

    ######################## Initialisation ############################################################

    Survie = Survieg = np.array([1])
    D = Dmoy = np.array([0])
    IncertSurvie = IncertSurvieg = Incert_dose = Edepmoy_n = np.array([])


    ######################## Root ##################################################################

    s, sg, dosen, dosec, dosem = (np.array([]) for i in range(5))
    nb_nucl_traversées_par_la_particule_tab_sur_toutes_simus, surviel_append_sur_une_simu,\
        surviel_append_sur_toutes_simus, survieg_append_sur_une_simu, survieg_append_sur_toutes_simus,\
        dosen_append_sur_une_simu, dosen_append_sur_toutes_simus, dosec_append_sur_une_simu, dosec_append_sur_toutes_simus,\
        dosem_append_sur_une_simu, dosem_append_sur_toutes_simus, dosen_c_append_sur_toutes_simus,\
        Ratio_CrossFire_Noyau_sur_toutes_simus, Ratio_CrossFire_Noyau_sur_toutes_simus_zone1,\
        Ratio_CrossFire_Noyau_sur_toutes_simus_zone2, Ei_Ef_sum_sur_toutes_simus,\
        Nombre_particules_par_noyau_sur_toutes_simus, TCP_append_sur_toutes_simus, TCP_test_formula_append_sur_toutes_simus,\
        SimulationId, Dose_Bio_append_sur_toutes_simus, test_spectre_diff =\
        ([] for nb_arrays in range(22))

    indexes_root_files_without_errors, indexes_root_files_without_errors_np, nb_files_with_errors = \
        id_deletion_of_root_outputs_with_errors()

    print("nb_cellules_reel : ", nb_cellules_reel)

    ##################################################

    for i in range(index_of_first_root_output, len(indexes_root_files_without_errors)):
        sum_dose_Noyau_CrossFire = sum_dose_Noyau_NonCrossFire = sum_dose_Noyau_NonCrossFire_zone1 = \
        sum_dose_Noyau_NonCrossFire_zone2 = sum_dose_Noyau_CrossFire_zone1 = sum_dose_Noyau_CrossFire_zone2 = 0

        hdepart = indexes_root_files_without_errors_np[i*nb_files_for_one_complete_simulation]
        hfin = hdepart + nb_files_for_one_complete_simulation -1

        dosen_append_sur_une_simu_np = np.zeros(nb_cellules_reel)
        dosec_append_sur_une_simu_np = np.zeros(nb_cellules_reel)
        dosem_append_sur_une_simu_np = np.zeros(nb_cellules_reel)

        n_unique_tot_sur_une_simu = np.zeros(nb_cellules_reel)
        Ei_Ef_unique_sur_une_simu = np.zeros(nb_cellules_reel)
        Nombre_particules_par_noyau = np.zeros(nb_cellules_reel)

        SimulationId.append(i) #working only if no multithreading

        print("Simu " + f"{i+1}" + " sur " + f"{(len(indexes_root_files_without_errors))}")

        for h in range((hdepart), (hfin+1), nb_group_of_cells_considered):
            root_file_name = []
            for ind_division_simus in range(0, nb_group_of_cells_considered):
                root_file_name.append("Root/" + "outputMultiCellulaire/" + dossier_root +\
                                    nom_fichier_root + f"{h + ind_division_simus}" + "_t0" + ".root")

            print("Root file name : ", root_file_name)

            data, ind_alphaplusplus, ind_alphaplus, ind_helium, data_alpha, ind_EndOfRun, data_EdepCell = \
                ([] for nb_arrays in range(7))
            indice_available_diffusion_info = indice_available_edep_sph_info = 0

            for ind_division_simus in range(0, nb_group_of_cells_considered):

                f = uproot.open(root_file_name[ind_division_simus])
                d1 = f['cell']['nameParticle'].array(library="np")
                d2 = f['cell']['Ei'].array(library="np")
                d3 = f['cell']['Ef'].array(library="np")
                d4 = f['cell']['ID_Cell'].array(library="np")
                d5 = f['cell']['Cellule_D_Emission'].array(library="np")
                d6 = f['cell']['eventID'].array(library="np")
                d7 = f['cell']['fEdepn'].array(library="np")
                d8 = f['cell']['fEdepc'].array(library="np")
                try:
                    d9 = f['cell']['indice_if_diffusion'].array(library="np")
                    indice_available_diffusion_info = 1
                except:
                    print("indice_if_diffusion not available on these data")
                    indice_available_diffusion_info = 0
                try:
                    d10 = f['cell']['fEdep_sph'].array(library="np")
                    indice_available_edep_sph_info = 1
                except:
                    print("fEdep_sph not available on these data")
                    indice_available_edep_sph_info = 0

                if (indice_available_diffusion_info == 0):
                    data.append(np.core.records.fromarrays([d1, d2, d3, d4, d5, d6, d7, d8],names='nameParticle, Ei, Ef, ID_Cell, Cellule_D_Emission, eventID, fEdepn, fEdepc'))
                elif (indice_available_diffusion_info == 1 and indice_available_edep_sph_info == 0):
                    data.append(np.core.records.fromarrays([d1, d2, d3, d4, d5, d6, d7, d8, d9],names='nameParticle, Ei, Ef, ID_Cell, Cellule_D_Emission, eventID, fEdepn, fEdepc, indice_if_diffusion'))
                elif (indice_available_diffusion_info == 1 and indice_available_edep_sph_info == 1):
                    data.append(np.core.records.fromarrays([d1, d2, d3, d4, d5, d6, d7, d8, d9, d10],names='nameParticle, Ei, Ef, ID_Cell, Cellule_D_Emission, eventID, fEdepn, fEdepc, indice_if_diffusion, fEdep_sph'))

                ind_alphaplusplus.append((data[ind_division_simus])["nameParticle"]=='alpha')
                ind_alphaplus.append((data[ind_division_simus])["nameParticle"]=='alpha+')
                ind_helium.append((data[ind_division_simus])["nameParticle"]=='helium')

                data_alpha.append((np.concatenate(((data[ind_division_simus])[ind_alphaplusplus[ind_division_simus]],(data[ind_division_simus])[ind_alphaplus[ind_division_simus]],(data[ind_division_simus])[ind_helium[ind_division_simus]]))))

                (data_alpha[ind_division_simus])["eventID"] += (nb_cellules_reel / nb_group_of_cells_considered)*ind_division_simus

                ind_EndOfRun.append((data[ind_division_simus])["nameParticle"] == 'EndOfRun')

                data_EdepCell.append((data[ind_division_simus])[ind_EndOfRun[ind_division_simus]])


            data_alpha=np.concatenate([row for row in data_alpha])

            ########################## Vérification diffusion aux bonnes énergies ###############################

            if (indice_available_diffusion_info == 1):

                unique_data_alpha_eventID = np.unique(data_alpha['eventID'], return_index=True)

                ind_diff_0 = ind_diff_1 = len_unique = 0

                indices_ab = unique_data_alpha_eventID[1]

                unique_data_alpha_ind_diff_corresponding_to_unique_event_id = np.take(data_alpha['indice_if_diffusion'], indices_ab)

                for i in range(0,len(unique_data_alpha_ind_diff_corresponding_to_unique_event_id)):
                    if (unique_data_alpha_ind_diff_corresponding_to_unique_event_id[i]==0):
                        ind_diff_0+=1
                        len_unique+=1
                    elif (unique_data_alpha_ind_diff_corresponding_to_unique_event_id[i] == 1):
                        ind_diff_1 += 1
                        len_unique += 1

                print("% d'event sans diffusion = ", ind_diff_0/len_unique)
                print("% d'event avec diffusion = ", ind_diff_1/len_unique)

            ####################################################################################################

            print()

            print(" data avec Ei < Ef : ", data_alpha[np.where(data_alpha["Ei"] < data_alpha["Ef"])])

            print()

            ####################### Modification des ID de CPOP ###################################

            ################ data_alpha #########################################

            print("len(data_alpha) = ",  len(data_alpha))

            for ind_modif_id in range(0,len(data_alpha)):
                index_ID_Cell = np.where(real_id_cells == data_alpha[ind_modif_id]["ID_Cell"])
                data_alpha[ind_modif_id]["ID_Cell"] = Perfect_ID_Cells[index_ID_Cell]

                index_Cellule_D_Emission = np.where(real_id_cells == data_alpha[ind_modif_id]["Cellule_D_Emission"])
                data_alpha[ind_modif_id]["Cellule_D_Emission"] = Perfect_ID_Cells[index_Cellule_D_Emission]


            ################ data_EdepCell ######################################

            if test_file_not_empty != 0:

                for ind_dose in range(0, nb_group_of_cells_considered):
                    elements_To_Remove = []
                    for ind_modif_id in range(0, len(data_EdepCell[ind_dose])):
                        if (((data_EdepCell[ind_dose])[ind_modif_id]["ID_Cell"]) in deleted_id_txt):
                            # elements_To_Remove.append((data_EdepCell[ind_dose])[ind_modif_id])
                            elements_To_Remove.append(ind_modif_id)
                    data_EdepCell[ind_dose] = np.delete(data_EdepCell[ind_dose], elements_To_Remove, 0)

            for ind_dose in range(0, nb_group_of_cells_considered):
                for ind_modif_id in range(0,len(data_EdepCell[ind_dose])):
                    index_ID_Cell = np.where(real_id_cells == (data_EdepCell[ind_dose])[ind_modif_id]["ID_Cell"])
                    (data_EdepCell[ind_dose])[ind_modif_id]["ID_Cell"] = Perfect_ID_Cells[index_ID_Cell]

            ####################### Histo d'énergie ########################################################

            Ei = data_alpha["Ei"] # Energy in keV
            Ef = data_alpha["Ef"]

            He = np.histogram(Ei, bins = bins, range = (1, Emax)) # [0] nbr of part in bins [1] energy of bins
            LET_He=np.zeros(len(He[0]))
            E_He=np.zeros(len(He[0]))
            for k in range(0,len(He[0])):
                E_He[k] = (He[1][k] + He[1][k + 1]) / 2.

            ################################################################################################

            n1 = number_of_lethal_events_for_alpha_traversals(dn1_dE_continous_pre_calculated)

            n = edep_n = edep_c = edep_m = 0

            n_tab = (n1(Ei) - n1(Ef))

            # print("############################# Calcul de dose #######################################")


            for ind_dose in range(0, nb_group_of_cells_considered):
                dosen_append_sur_une_simu_np += (((data_EdepCell[ind_dose])["fEdepn"]) * KEV_IN_J / masses_nuclei)
                dosec_append_sur_une_simu_np += (((data_EdepCell[ind_dose])["fEdepc"]) * KEV_IN_J / masses_cytoplasms)


            count_cell_id = np.bincount((data_alpha["ID_Cell"]).astype(int))
            Ei_Ef_unique = np.bincount(data_alpha["ID_Cell"].astype(int), weights= Ei - Ef)

            while (len(Ei_Ef_unique) < (nb_cellules_reel)):
                Ei_Ef_unique = np.append(Ei_Ef_unique, 0)

            while (len(count_cell_id) < (nb_cellules_reel)):
                count_cell_id = np.append(count_cell_id, 0)

            Ei_Ef_unique_sur_une_simu += Ei_Ef_unique
            Nombre_particules_par_noyau += count_cell_id


            #################################### Cross-fire dose au noyau ##############################################################


            ############################################

            ind_NonCrossFire = data_alpha["ID_Cell"] == data_alpha["Cellule_D_Emission"]
            ind_CrossFire = data_alpha["ID_Cell"] != data_alpha["Cellule_D_Emission"]

            if (indice_available_diffusion_info == 1):
                ind_NonCrossFire = ((data_alpha["ID_Cell"] == data_alpha["Cellule_D_Emission"]) & (data_alpha["indice_if_diffusion"]==0))
                ind_CrossFire = ((data_alpha["ID_Cell"] != data_alpha["Cellule_D_Emission"]) & (data_alpha["indice_if_diffusion"]==1))

            data_Noyau_NonCrossFire = data_alpha[ind_NonCrossFire]

            dose_noyau_NonCrossFire=data_Noyau_NonCrossFire["Ei"]-data_Noyau_NonCrossFire["Ef"]

            dose_Noyau_CrossFire=data_alpha["Ei"]-data_alpha["Ef"]

            dose_Noyau_CrossFire=np.setdiff1d(dose_Noyau_CrossFire,dose_noyau_NonCrossFire)

            sum_dose_Noyau_CrossFire+=np.sum(dose_Noyau_CrossFire)
            sum_dose_Noyau_NonCrossFire+=np.sum(dose_noyau_NonCrossFire)

            ################################################
            indice_zone1 = zone_cell == 1
            indice_zone2 = zone_cell == 2

            Ei_Ef_unique_NonCrossFire = np.bincount(((data_alpha["ID_Cell"])[ind_NonCrossFire]).astype(int), weights=Ei[ind_NonCrossFire] - Ef[ind_NonCrossFire])
            Ei_Ef_unique_CrossFire = np.bincount(((data_alpha["ID_Cell"])[ind_CrossFire]).astype(int), weights=Ei[ind_CrossFire] - Ef[ind_CrossFire])

            while (len(Ei_Ef_unique_NonCrossFire) < (nb_cellules_reel)):
                Ei_Ef_unique_NonCrossFire = np.append(Ei_Ef_unique_NonCrossFire, 0)

            while (len(Ei_Ef_unique_CrossFire) < (nb_cellules_reel)):
                Ei_Ef_unique_CrossFire = np.append(Ei_Ef_unique_CrossFire, 0)

            Ei_Ef_unique_NonCrossFire_zone1 = Ei_Ef_unique_NonCrossFire[indice_zone1]
            Ei_Ef_unique_NonCrossFire_zone2 = Ei_Ef_unique_NonCrossFire[indice_zone2]

            Ei_Ef_unique_CrossFire_zone1 = Ei_Ef_unique_CrossFire[indice_zone1]
            Ei_Ef_unique_CrossFire_zone2 = Ei_Ef_unique_CrossFire[indice_zone2]

            sum_dose_Noyau_NonCrossFire_zone1 += np.sum(Ei_Ef_unique_NonCrossFire_zone1)
            sum_dose_Noyau_NonCrossFire_zone2 += np.sum(Ei_Ef_unique_NonCrossFire_zone2)

            sum_dose_Noyau_CrossFire_zone1 += np.sum(Ei_Ef_unique_CrossFire_zone1)
            sum_dose_Noyau_CrossFire_zone2 += np.sum(Ei_Ef_unique_CrossFire_zone2)


            #################################### Nombre de cellules traversées par particule ###########################################

            nb_particules_tot=len(np.unique(data_alpha["eventID"]))

            count_event_id = np.bincount((data_alpha["eventID"]).astype(int))

            for id_part in range(0, len(np.unique(data_alpha["eventID"]))):
                nb_nucl_traversées_par_la_particule=count_event_id[id_part]
                nb_nucl_traversées_par_la_particule_tab_sur_toutes_simus.append(nb_nucl_traversées_par_la_particule)


            # print("############################# Calcul de survie #######################################")

            n_unique=np.bincount(data_alpha["ID_Cell"].astype(int), weights=n_tab)
            while (len(n_unique) < (nb_cellules_reel)):
                n_unique=np.append(n_unique,0)

            n_unique_tot_sur_une_simu+=n_unique

        print()

        print(" data avec Ei < Ef : ", data_alpha[np.where(data_alpha["Ei"]<data_alpha["Ef"])])

        print()

        sum_dose_Noyau_tot=sum_dose_Noyau_CrossFire+sum_dose_Noyau_NonCrossFire
        Ratio_CrossFire_Noyau_sur_une_simu=sum_dose_Noyau_CrossFire/sum_dose_Noyau_tot
        Ratio_CrossFire_Noyau_sur_toutes_simus.append(Ratio_CrossFire_Noyau_sur_une_simu)

        Ratio_CrossFire_Noyau_sur_une_simu_zone1=0
        Ratio_CrossFire_Noyau_sur_toutes_simus_zone1.append(Ratio_CrossFire_Noyau_sur_une_simu_zone1)

        Ratio_CrossFire_Noyau_sur_une_simu_zone2=0
        Ratio_CrossFire_Noyau_sur_toutes_simus_zone2.append(Ratio_CrossFire_Noyau_sur_une_simu_zone2)

        Ei_Ef_sum_sur_toutes_simus.append(Ei_Ef_unique_sur_une_simu)
        Nombre_particules_par_noyau_sur_toutes_simus.append(Nombre_particules_par_noyau)

        surviel_append_sur_une_simu=np.exp(-n_unique_tot_sur_une_simu)

        surviel_append_sur_une_simu[np.where(surviel_append_sur_une_simu == 0)] = 10**(-299)

        surviel_append_sur_toutes_simus.append(surviel_append_sur_une_simu)

        alpha_ref = 0.313  # HSG
        beta_ref = 0.0615  # HSG

        print()

        Dose_Bio_append_sur_une_simu = (np.sqrt(alpha_ref ** 2 - 4 * beta_ref * np.log(surviel_append_sur_une_simu)) - alpha_ref) / (2 * beta_ref)
        Dose_Bio_append_sur_toutes_simus.append(Dose_Bio_append_sur_une_simu)

        exp_surviel = np.exp(-np.asarray(surviel_append_sur_une_simu))
        TCP_une_simu = np.prod(exp_surviel)
        TCP_test_formula = np.prod(1 - surviel_append_sur_une_simu)
        TCP_append_sur_toutes_simus.append(TCP_une_simu)
        TCP_test_formula_append_sur_toutes_simus.append(TCP_test_formula)

        survieg_append_sur_une_simu=np.exp(-n_unique_tot_sur_une_simu - BETAG[type_cell] * (dosen_append_sur_une_simu_np ** 2))
        survieg_append_sur_toutes_simus.append(survieg_append_sur_une_simu)


        dosen_append_sur_toutes_simus.append(dosen_append_sur_une_simu_np)
        dosec_append_sur_toutes_simus.append(dosec_append_sur_une_simu_np)
        dosem_append_sur_toutes_simus.append(dosem_append_sur_une_simu_np)
        dosen_c_append_sur_toutes_simus.append(dosen_append_sur_une_simu_np+dosec_append_sur_une_simu_np)

        Dose_Spheroid = data_EdepCell[0]["fEdep_sph"]* KEV_IN_J / masse_tum

        ############################################################################################################################

    nb_nucl_traversées_par_la_particule_tab_sur_toutes_simus_np=np.asarray(nb_nucl_traversées_par_la_particule_tab_sur_toutes_simus)
    surviel_append_sur_toutes_simus_np=np.asarray(surviel_append_sur_toutes_simus)
    survieg_append_sur_toutes_simus_np=np.asarray(survieg_append_sur_toutes_simus)


    Ei_Ef_sum_sur_toutes_simus_np = np.asarray(Ei_Ef_sum_sur_toutes_simus)
    Nombre_particules_par_noyau_sur_toutes_simus_np = np.asarray(Nombre_particules_par_noyau_sur_toutes_simus)
    Nombre_particules_par_noyau_sur_toutes_simus_sans_zero_np = np.array([])
    Ei_Ef_sum_sur_toutes_simus_sans_zero_np = np.array([])

    ########### Suppression des cells qui n'ont pas été traversées par des alpha, pour le calcul d'Edep moy quand les particules touchent des noyaux ###########

    nb_cells_suppr = 0

    for i in range(0,(nb_complete_simulations - nb_files_with_errors)):
        ind_Noyau_Avec_Zero_Particules = [index for index, value in enumerate(Nombre_particules_par_noyau_sur_toutes_simus_np[i]) if value == 0]

        Nombre_particules_par_noyau_sur_toutes_simus_sans_zero_np = np.append(Nombre_particules_par_noyau_sur_toutes_simus_sans_zero_np, np.delete(Nombre_particules_par_noyau_sur_toutes_simus_np[i], ind_Noyau_Avec_Zero_Particules))
        Ei_Ef_sum_sur_toutes_simus_sans_zero_np = np.append(Ei_Ef_sum_sur_toutes_simus_sans_zero_np, np.delete(Ei_Ef_sum_sur_toutes_simus_np[i], ind_Noyau_Avec_Zero_Particules))

        nb_cells_suppr += len(ind_Noyau_Avec_Zero_Particules)


    print()

    Edep_dans_noy_par_particule_np = Ei_Ef_sum_sur_toutes_simus_sans_zero_np/Nombre_particules_par_noyau_sur_toutes_simus_sans_zero_np
    mean_Edep_dans_noy_par_particule = np.mean(Edep_dans_noy_par_particule_np)

    dosen_append_sur_toutes_simus_np=np.asarray(dosen_append_sur_toutes_simus)
    dosec_append_sur_toutes_simus_np=np.asarray(dosec_append_sur_toutes_simus)
    dosem_append_sur_toutes_simus_np=np.asarray(dosem_append_sur_toutes_simus)
    dosen_c_append_sur_toutes_simus_np=np.asarray(dosen_c_append_sur_toutes_simus)
    Ratio_CrossFire_Noyau_sur_toutes_simus_np=np.asarray(Ratio_CrossFire_Noyau_sur_toutes_simus)
    Ratio_CrossFire_Noyau_sur_toutes_simus_zone1_np=np.asarray(Ratio_CrossFire_Noyau_sur_toutes_simus_zone1)
    Ratio_CrossFire_Noyau_sur_toutes_simus_zone2_np=np.asarray(Ratio_CrossFire_Noyau_sur_toutes_simus_zone2)
    TCP_append_sur_toutes_simus_np = np.asarray(TCP_append_sur_toutes_simus)
    TCP_test_formula_append_sur_toutes_simus_np = np.asarray(TCP_test_formula_append_sur_toutes_simus)
    Dose_Bio_append_sur_toutes_simus_np = np.asarray(Dose_Bio_append_sur_toutes_simus)

    surviel=np.mean(surviel_append_sur_toutes_simus_np,axis=0)
    survieg=np.mean(survieg_append_sur_toutes_simus_np,axis=0)
    dosen_tot=np.mean(dosen_append_sur_toutes_simus_np,axis=0)
    dosec_tot=np.mean(dosec_append_sur_toutes_simus_np,axis=0)
    dosem_tot=np.mean(dosem_append_sur_toutes_simus_np,axis=0)
    dosen_c_tot=np.mean(dosen_c_append_sur_toutes_simus_np,axis=0)
    TCP = np.mean(TCP_append_sur_toutes_simus_np, axis=0)
    TCP_formula2 = np.mean(TCP_test_formula_append_sur_toutes_simus_np,axis=0)
    dose_bio = np.mean(Dose_Bio_append_sur_toutes_simus_np, axis=0)
    mean_CrossFire= np.full(nb_cellules_reel ,np.mean(Ratio_CrossFire_Noyau_sur_toutes_simus_np)*100)

    dose_tot=dosen_tot+dosec_tot

    sum_dose_noyau = np.mean(np.sum(dosen_append_sur_toutes_simus_np,axis=1))
    sum_dose_tot = np.mean(np.sum(dosen_c_append_sur_toutes_simus_np,axis=1))

    incertdmoy_n_tot=np.std(dosen_append_sur_toutes_simus_np,axis=0)
    incertdmoy_c_tot=np.std(dosec_append_sur_toutes_simus_np,axis=0)
    incertdmoy_m_tot=np.std(dosem_append_sur_toutes_simus_np,axis=0)
    incertdmoy_n_c_tot=np.std(dosen_c_append_sur_toutes_simus_np,axis=0)
    incert_crossfire=np.full(nb_cellules_reel ,np.std(Ratio_CrossFire_Noyau_sur_toutes_simus_np)*100)
    incert_mean_cell_survival=np.full(nb_cellules_reel ,np.std(np.mean(surviel_append_sur_toutes_simus_np,axis=1)))
    incert_dose_bio=np.std(np.mean(Dose_Bio_append_sur_toutes_simus_np, axis=1))

    print()
    print("nb de config avec erreur", nb_files_with_errors)
    print()

    alpha=0.5

    SF=np.sum(np.exp(-alpha*dosen_c_tot))/nb_cellules_reel

    EUD=-np.log(SF)/alpha

    # print(TCP_append_sur_toutes_simus_np)

    print("TCP =")
    print(TCP, "+-", 2*np.std(TCP_append_sur_toutes_simus_np)/np.sqrt(nb_complete_simulations))

    print("TCP test formula = ", TCP_formula2 , "+-", 2*np.std(TCP_test_formula_append_sur_toutes_simus_np)/np.sqrt(nb_complete_simulations))

    print("EUD =")
    print(EUD)

    print("Min dose totale absorbée=")
    print(np.mean(np.min(dosen_c_append_sur_toutes_simus_np,axis=1)), "+-", 2*np.std(np.min(dosen_c_append_sur_toutes_simus_np,axis=1)))

    print("Max dose totale absorbée=")
    print(np.mean(np.max(dosen_c_append_sur_toutes_simus_np,axis=1)), "+-", 2*np.std(np.max(dosen_c_append_sur_toutes_simus_np,axis=1)))

    # Les 2 méthodes de calcul de la moyenne sont les mêmes, mais la seconde calcule correctement l'incertitude

    print("Moyenne des doses absorbées aux noyaux=")
    print(np.mean(np.mean(dosen_append_sur_toutes_simus_np,axis=1)), "+-", 2*np.std(np.mean(dosen_append_sur_toutes_simus_np,axis=1)))

    print("Moyenne des doses absorbées aux cytoplasmes=")
    print(np.mean(np.mean(dosec_append_sur_toutes_simus_np,axis=1)), "+-", 2*np.std(np.mean(dosec_append_sur_toutes_simus_np,axis=1)))

    print("Moyenne des doses absorbées aux cellules=")
    print(np.mean(np.mean(dosen_c_append_sur_toutes_simus_np,axis=1)), "+-", 2*np.std(np.mean(dosen_c_append_sur_toutes_simus_np,axis=1)))

    print("Sum dose cellules = ")
    print(np.sum(dosen_c_tot))

    print("Energie moyenne déposée par une particule quand elle touche un noyau=")
    print(mean_Edep_dans_noy_par_particule, "+-", np.std(Edep_dans_noy_par_particule_np)/np.sqrt(len(Edep_dans_noy_par_particule_np)))

    print("Nombre moyen de noyaux traversés par une particule en moyenne=")
    print(np.mean(nb_nucl_traversées_par_la_particule_tab_sur_toutes_simus), "+-", np.std(nb_nucl_traversées_par_la_particule_tab_sur_toutes_simus)/np.sqrt(len(nb_nucl_traversées_par_la_particule_tab_sur_toutes_simus))) #Pour distribution, voir histo

    print("Cross-fire dose au noyau en moyenne =")
    print(np.mean(Ratio_CrossFire_Noyau_sur_toutes_simus_np)*100, "%", "+-", 2*np.std(Ratio_CrossFire_Noyau_sur_toutes_simus_np)*100, "%")

    print("Cross-fire dose au noyau en moyenne, dans la zone 1 =")
    print(np.mean(Ratio_CrossFire_Noyau_sur_toutes_simus_zone1_np)*100, "%", "+-", 2*np.std(Ratio_CrossFire_Noyau_sur_toutes_simus_zone1_np)*100, "%")

    print("Cross-fire dose au noyau en moyenne, dans la zone 2 =")
    print(np.mean(Ratio_CrossFire_Noyau_sur_toutes_simus_zone2_np)*100, "%", "+-", 2*np.std(Ratio_CrossFire_Noyau_sur_toutes_simus_zone2_np)*100, "%")

    print("Moyenne des survies cellulaires=")
    print(np.mean(np.mean(surviel_append_sur_toutes_simus_np,axis=1)), "+-", 2*np.std(np.mean(surviel_append_sur_toutes_simus_np,axis=1)))

    print("Moyenne des doses biologiques")
    print(np.mean(np.mean(Dose_Bio_append_sur_toutes_simus_np, axis=1)), "+-", 2 * np.std(np.mean(Dose_Bio_append_sur_toutes_simus_np, axis=1)))

    print("Max des doses biologiques")
    print(np.mean(np.max(Dose_Bio_append_sur_toutes_simus_np, axis=1)), "+-",2 * np.std(np.max(Dose_Bio_append_sur_toutes_simus_np, axis=1)))

    if(indice_available_edep_sph_info):
        print("Dose sphéroïde")
        print(Dose_Spheroid[0], "Gy")

    print("Nombre de simus fonctionnelles = ")
    print(nb_complete_simulations - nb_files_with_errors)

    print()

    id_cell_arr=np.arange(nb_cellules_reel)


    fname = "AnalysisResults/" + study_type_folder_name + "/" + nom_dossier_pour_excel_analyse + "/" + "Emission" + cell_compartment_radionuclide + ".xlsx"
    name_columns = ['ID_Cell','Zone_Cell' ,'Survie locale', 'Incert Survie Locale','Survie globale', 'Dmoy au noyau, par simu (Gy)','Incert Dmoy au noyau par simu(Gy)','Dmoy cyto, par simu (Gy)','incert Dmoy cyto par simu (Gy)','Dmoy au noyau+cyto, par simu (Gy)','incert Dmoy noy+cyto par simu (Gy)', 'Ratio Crossfire %','Incert Ratio Crossfire %','Ratio Crossfire zone 1 %', 'Ratio Crossfire zone 2 %','Nb noyau par particule', 'Sum dose sur tous noyaux', 'Sum dose sur toutes cellules', 'Edep moy par particule dans noyau', 'Incert Edep moy par particule dans noyau', 'Nb noyaux par particule', 'Incert nb noyaux par particule' ,'TCP', 'Incert_TCP','TCP formula 2','Incert TCP formula 2','EUD','Dose bio','IncertDoseBio', "DoseSpheroid"]
    results_to_write = np.array(np.concatenate([[name_columns],np.transpose([id_cell_arr,zone_cell,surviel, incert_mean_cell_survival, survieg,dosen_tot,np.full(nb_cellules_reel,incertdmoy_n_tot),dosec_tot,np.full(nb_cellules_reel,incertdmoy_c_tot),dosen_c_tot, np.full(nb_cellules_reel,incertdmoy_n_c_tot), mean_CrossFire, incert_crossfire,np.full(nb_cellules_reel,0), np.full(nb_cellules_reel,0), np.full(nb_cellules_reel, np.full(nb_cellules_reel, np.mean(nb_nucl_traversées_par_la_particule_tab_sur_toutes_simus))), np.full(nb_cellules_reel,sum_dose_noyau), np.full(nb_cellules_reel,sum_dose_tot), np.full(nb_cellules_reel,mean_Edep_dans_noy_par_particule), np.full(nb_cellules_reel, np.std(Edep_dans_noy_par_particule_np)/np.sqrt(len(Edep_dans_noy_par_particule_np))), np.full(nb_cellules_reel, np.mean(nb_nucl_traversées_par_la_particule_tab_sur_toutes_simus)), np.full(nb_cellules_reel, np.std(nb_nucl_traversées_par_la_particule_tab_sur_toutes_simus)/np.sqrt(len(nb_nucl_traversées_par_la_particule_tab_sur_toutes_simus))),np.full(nb_cellules_reel,TCP), np.full(nb_cellules_reel, np.std(TCP_append_sur_toutes_simus_np)),np.full(nb_cellules_reel,TCP_formula2), np.full(nb_cellules_reel, np.std(TCP_test_formula_append_sur_toutes_simus_np)),np.full(nb_cellules_reel,EUD), dose_bio, np.full(nb_cellules_reel,incert_dose_bio), np.full(nb_cellules_reel,Dose_Spheroid[0])])]))
    wb = Workbook()
    wb.new_sheet("AnalysisResults", data=results_to_write)
    wb.save(fname)


    fname = "AnalysisResults/" + study_type_folder_name + "/" + nom_dossier_pour_excel_analyse + "/" + "Emission" + cell_compartment_radionuclide + "SimID" + ".xlsx"
    name_columns = ['SimulationID', 'TCP', 'Dose moyenne noyau', 'Incert dose moyenne noyau', 'Dose moyenne cellule', 'Incert dose moyenne cellule', 'Dose bio moyenne', 'Incert dose bio moyenne']
    results_to_write = np.array(np.concatenate([[name_columns],np.transpose([SimulationId, TCP_append_sur_toutes_simus_np, np.mean(dosen_append_sur_toutes_simus_np,axis=1), np.std(dosen_append_sur_toutes_simus_np,axis=1)/np.sqrt(nb_cellules_reel), np.mean(dosen_c_append_sur_toutes_simus_np,axis=1), np.std(dosen_c_append_sur_toutes_simus_np,axis=1)/np.sqrt(nb_cellules_reel), np.mean(Dose_Bio_append_sur_toutes_simus_np,axis=1), np.std(Dose_Bio_append_sur_toutes_simus_np,axis=1)/np.sqrt(nb_cellules_reel)])]))
    wb = Workbook()
    wb.new_sheet("AnalysisResults_SimID", data=results_to_write)
    wb.save(fname)

    print(" Temps total =  ", (time.time() - START_TIME)//60, "minutes", (time.time() - START_TIME)%60, "secondes")

def add_new_buttons_to_graphic_window():
    global r_sph, nom_config, spheroid_compaction, xml_geom, nb_cellules_xml, cell_compartment_radionuclide,\
    nb_complete_simulations, nb_files_for_one_complete_simulation, nb_group_of_cells_considered, SimulationName, study_type_folder_name,\
    bool_diff, rn_name, nb_particles_per_cell, type_cell, available_data_date, available_data_name_file, \
    data_cb, nom_fichier_root

    geom_list = ["Elg030um75CP", "Elg050um75CP", "Elg070um75CP", "Elg160um75CP", "Elg095um25CP",
                 "Elg095um50CP", "Elg095um75CP", "Elg095um75CP_2", "Elg100um40CP"]
    cp_list= [25, 50, 75]
    r_sph = geom_list[geom_cb.current()][3:6]
    print("r_sph : ", r_sph)
    print(float(r_sph))

    nom_config = (geom_list[geom_cb.current()])  # Les fichiers contenant les masses de toutes les cellules, et ceux des ID de cellules supprimés de CPOP à G4, sont appelés MassesCell_nom_config.txt, et IDCell_nom_config.txt

    # cp = (cp_list[geom_cb.current()])
    spheroid_compaction = geom_list[geom_cb.current()][8:10]
    print("spheroid_compaction : ", spheroid_compaction)

    xml_geom = "Cpop_Geom_XML/" + nom_config + ".cfg" + ".xml"

    nb_cellules_xml = geometry_informations.count_number_of_cells_in_xml_file(xml_geom)  # Nombre de cellules contenues dans le fichier .xml de géométrie créé par CPOP
    print("nb_cellules_xml", nb_cellules_xml)

    cell_compartment_radionuclide = (distrib_name_cb.get())

    study_type = radiovalue_study_type.get()  # 0 for internalization study, 1 for labeling study

    nb_complete_simulations = int(nb_simulations_entry.get())
    nb_files_for_one_complete_simulation = 1    # Corresponds to the number of jobs/root outputs in which one complete
                                                # simulation is divided. Usually one but can be changed for very long
                                                # simulations.
    nb_group_of_cells_considered = 1  #Usually 1 but can be changed for very long simulations.
                                      # E.g. : if = 2, half the jobs were sent in 50% of the cells
                                      #and the other jobs in the other 50% of cells.

    if (study_type == 0):
        SimulationName = cell_compartment_radionuclide
        study_type_folder_name = "Internalization"
    elif (study_type == 1):
        LabelingPercentage_get = Labeling_Percentage_Entry.get()
        LabelingPercentage_name = str(LabelingPercentage_get) + "_Percent"
        SimulationName = LabelingPercentage_name
        study_type_folder_name = "Labeling"

    # output_folders = glob.glob("Root/outputMultiCellulaire/" + study_type_folder_name + "/", recursive = True)

    output_path = "Root/outputMultiCellulaire/" + study_type_folder_name + "/"

    output_folders_name = [f for f in os.listdir(output_path)]

    print(output_folders_name)

    bool_diff = ["Yes","No"]
    rn_name = radionuclide_name_entry.get()
    nb_particles_per_cell = ["1", "2", "3", "4", "5" ,"6", "7", "8", "9" ,"10", "42"]

    type_cell = cell_line_cb.current()
    print("type cell is number : ", type_cell)

    available_data_date = []
    available_data_name_file = []
    for i in range(0, len(output_folders_name)):
        if ("_" + cell_compartment_radionuclide + "_" + str(spheroid_compaction) + "CP_" + str(r_sph) + "um_" + rn_name + "_diff" + bool_diff[diffusion_list.current()] + "_" + str(nb_particles_per_cell[number_particles_per_cell_list.current()] + "ppc")) in \
                output_folders_name[i]:
            # available_data.append(output_folders_name[i])
            available_data_date.append(output_folders_name[i][0:10])
            available_data_name_file.append(output_folders_name[i])

    print(available_data_name_file)

    # window_data = tkinter.Toplevel()
    # window_data.geometry("700x300")

    List_data_Txt = tkinter.Label(window, text="Data available :", fg='blue')
    List_data_Txt.place(x=100, y=510)
    selected_data = tkinter.StringVar()
    data_cb = tkinter.ttk.Combobox(window, width=35, textvariable=selected_data)
    data_cb['values'] = available_data_date
    data_cb.place(x=400, y=510)

    nom_fichier_root = "output_"  # Les fichiers root, contenus dans le dossier_root, s'appellent nom_fichier_root{0,...}.root

    Validate_data = tkinter.Button(window, text="Validate", command=main)
    Validate_data.place(x=480, y=560)

graphic_window()
