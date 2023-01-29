#!/usr/bin/env python3
"""
Script allowing to convert .root raw data of Geant4 in data of interest

Returns
-------
    doses to cell nucleus and cytoplams
    cell survivals
    cross-fire information
"""

import geometry_informations
import math
import numpy as np
import os
import pandas as pd
from pyexcelerate import Workbook
import scipy.integrate
import scipy.interpolate as interpolate
import sys
import time
import tkinter
import tkinter.ttk
import uproot

KEV_IN_J = 1.60218 * 1e-16
WATER_DENSITY = 1e-15  #kg/µm³
UNIT_COEFFICIENT_A = KEV_IN_J / WATER_DENSITY # Gy.µm³.keV-1

SIG0 = [49*np.pi, 24.01*np.pi, 34.81*np.pi] #From Mario calculations
energies_valid_for_alpha_beta_approximation = np.arange(200,90001)

################ Fit parameters used in article TCP RIV-alpha ###################
# Y0 = [0.06072969, 0.02562553, 0.03934994] #HSG, V79, CHO-K1
# A = [-0.18385472, -0.10426184, -0.11163773]
# W = [3.05093045, 2.87758559, 3.20398251]
# XC = [0.46545609, 0.38084839, 0.48452192]

################ Fit parameters from 2022/12/16 ################# TO DO : evaluate differences with previous parameters
Y0 = [0.06486164, 0.02722410, 0.04221387] #HSG, V79, CHO-K1
A = [-0.26336407, -0.11801719, -0.19357751]
W = [3.39940424, 2.97713123, 3.90866411]
XC = [-0.00863166, 0.23348883, -0.25238105]

BETAG = [0.0961, 0.0405, 0.0625]  # constante de Monini et al. 2019


Emax=8000 #Energie max des ions Hélium émis, en keV

radius_cell_line = 7  # Rayon du noyau de la lignée HSG, en µm
surface_centerslice_cell_line = math.pi * radius_cell_line ** 2 #µm²
length_of_cylinderslice_cell = 1 #µm

bins = 200
start_time = time.perf_counter()

np.set_printoptions(threshold=sys.maxsize)
#np.set_printoptions(threshold = False)

global labeling_percentage_entry, cell_compartment_label, cell_compartment_combobox, verbose

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

    tables_conversion_energy_in_let = pandas.read_excel(f"E_TEL/conversion_tables_{data_base}.xlsx").to_records()
    energy_list = tables_conversion_energy_in_let['E(keV)']
    corresponding_let_list = tables_conversion_energy_in_let['LET(keV/um)']
    continuous_function_to_convert_energy_in_let = interpolate.interp1d(energy_list, corresponding_let_list,
                                                                        fill_value="extrapolate", kind= "linear")

    return continuous_function_to_convert_energy_in_let(energy)

def beta_nanox(energy,type_cell): #energy in MeV/nucleon
    return Y0[type_cell] + (A[type_cell]/(2*np.pi))*(W[type_cell]/(((energy-XC[type_cell])**2)+(W[type_cell]**2)/4))

def alpha_nanox(energy,type_cell):
    conv_let_e_srim = conversion_energy_in_let("SRIM", energy)
    b = BETAG[type_cell] * (((UNIT_COEFFICIENT_A * conv_let_e_srim * 0.8)/SIG0[type_cell]) ** 2)
    return (SIG0[type_cell] +
           ((UNIT_COEFFICIENT_A * conv_let_e_srim * (b-1))*np.sqrt(beta_nanox(energy/4000,type_cell)/(b+(b*b)/2))))\
           /(UNIT_COEFFICIENT_A * conv_let_e_srim)

def dn1_de_continous():
    conversion_energy_in_let_srim_alpha_beta_approximation_range = \
        conversion_energy_in_let("SRIM", energies_valid_for_alpha_beta_approximation)
    conversion_energy_in_let_g4_alpha_beta_approximation_range = \
        conversion_energy_in_let("G4", energies_valid_for_alpha_beta_approximation)
    dn1_de = -np.log(1 - alpha_nanox(energies_valid_for_alpha_beta_approximation,type_cell) \
                     * UNIT_COEFFICIENT_A*conversion_energy_in_let_srim_alpha_beta_approximation_range \
                     / surface_centerslice_cell_line) \
                     / (length_of_cylinderslice_cell * conversion_energy_in_let_g4_alpha_beta_approximation_range)
                    #calculation of number of lethal events per keV, via Mario's approximations
    dn1_de_interpolated = interpolate.interp1d(energies_valid_for_alpha_beta_approximation,
                                               dn1_de, fill_value="extrapolate", kind="linear")
    return dn1_de_interpolated

def number_of_lethal_events_for_alpha_traversals(dn1_de_function):
    """
    Returns the function that converts an energy E into the cumulated number of lethal damage from 0 to E
    """
    energie_table_binned = np.linspace(0, Emax, num=bins)
    f_he_cumulative_int = scipy.integrate.cumtrapz(dn1_de_function(energie_table_binned),
                                                   energie_table_binned, initial=0)
    n1 = interpolate.interp1d(energie_table_binned, f_he_cumulative_int, fill_value="extrapolate",
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
        if positions_cell[index] < radius_zone_1:
            zone_cell[index] = 1
            nb_cell_zone_1 += 1
        elif positions_cell[index] < radius_zone_2:
            zone_cell[index] = 2
            nb_cell_zone_2 += 1
    return zone_cell, nb_cell_zone_1, nb_cell_zone_2

def calculations_from_root_file(analysis_results, simulation_id):
    """
    Opens root file corresponding to a MC simulation and calculates quantities like cell survivals
    Returns Check
    """

    analysis_dataframe = pd.DataFrame()
    analysis_dataframe['simulation_id'] = simulation_id

    root_file_name = f"Root/outputMultiCellulaire/{dossier_root}{nom_fichier_root}{i}_t0.root"

    if verbose == 1:
        print("Root file name : ", root_file_name)

    ind_alphaplusplus, ind_alphaplus, ind_helium, data_alpha, ind_end_of_run, data_edep_cell = \
        ([] for _ in range(6))
    indice_available_diffusion_info = indice_available_edep_sph_info = 0

    f = uproot.open(root_file_name)
    name_particle_root = f['cell']['nameParticle'].array(library="np")

    ei_root = f['cell']['Ei'].array(library="np")
    ef_root = f['cell']['Ef'].array(library="np")
    id_cell_root = f['cell']['ID_Cell'].array(library="np")
    emission_cell_root = f['cell']['Cellule_D_Emission'].array(library="np")
    event_id_root = f['cell']['eventID'].array(library="np")
    energy_deposited_nucleus_root = f['cell']['fEdepn'].array(library="np")
    energy_deposited_cytoplasm_root = f['cell']['fEdepc'].array(library="np")
    try:
        diffusion_index_root = f['cell']['indice_if_diffusion'].array(library="np")
        indice_available_diffusion_info = 1
    except:
        print("indice_if_diffusion not available on these data")
        indice_available_diffusion_info = 0
    try:
        energy_deposited_spheroid_root = f['cell']['fEdep_sph'].array(library="np")
        indice_available_edep_sph_info = 1
    except:
        print("fEdep_sph not available on these data")
        indice_available_edep_sph_info = 0

    if indice_available_diffusion_info == 0:
        root_data_opened = np.core.records.fromarrays([name_particle_root, ei_root, ef_root,
                                                            id_cell_root, emission_cell_root,
                                                            event_id_root, energy_deposited_nucleus_root,
                                                            energy_deposited_cytoplasm_root],
                                                           names='nameParticle, Ei, Ef, ID_Cell, Cellule_D_Emission, eventID, fEdepn, fEdepc')
    elif indice_available_diffusion_info == 1 and indice_available_edep_sph_info == 0:
        root_data_opened = np.core.records.fromarrays([name_particle_root, ei_root, ef_root,
                                                            id_cell_root, emission_cell_root,
                                                            event_id_root, energy_deposited_nucleus_root,
                                                            energy_deposited_cytoplasm_root,
                                                            diffusion_index_root],
                                                           names='nameParticle, Ei, Ef, ID_Cell, Cellule_D_Emission, eventID, fEdepn, fEdepc,'
                                                                 ' indice_if_diffusion')
    elif indice_available_diffusion_info == 1 and indice_available_edep_sph_info == 1:
        root_data_opened = np.core.records.fromarrays([name_particle_root, ei_root, ef_root,
                                                            id_cell_root, emission_cell_root,
                                                            event_id_root, energy_deposited_nucleus_root,
                                                            energy_deposited_cytoplasm_root,
                                                            diffusion_index_root,
                                                            energy_deposited_spheroid_root],
                                                           names='nameParticle, Ei, Ef, ID_Cell, Cellule_D_Emission, eventID, fEdepn, fEdepc,'
                                                                 ' indice_if_diffusion, fEdep_sph')

    ind_alphaplusplus = root_data_opened["nameParticle"] == 'alpha'
    ind_alphaplus = root_data_opened["nameParticle"] == 'alpha+'
    ind_helium = root_data_opened["nameParticle"] == 'helium'

    data_alpha = (np.concatenate((root_data_opened[ind_alphaplusplus],
                                       root_data_opened[ind_alphaplus],
                                       root_data_opened[ind_helium])))

    ind_end_of_run = root_data_opened["nameParticle"] == 'EndOfRun'

    data_edep_cell = root_data_opened[ind_end_of_run]

    # data_alpha = np.concatenate([row for row in data_alpha])

    ########################## Vérification diffusion aux bonnes énergies ###############################

    if indice_available_diffusion_info == 1:

        unique_data_alpha_event_id = np.unique(data_alpha['eventID'], return_index=True)

        ind_diff_0 = ind_diff_1 = len_unique = 0

        indices_ab = unique_data_alpha_event_id[1]

        unique_data_alpha_ind_diff_corresponding_to_unique_event_id = \
            np.take(data_alpha['indice_if_diffusion'], indices_ab)

        for i in range(0, len(unique_data_alpha_ind_diff_corresponding_to_unique_event_id)):
            if unique_data_alpha_ind_diff_corresponding_to_unique_event_id[i] == 0:
                ind_diff_0 += 1
                len_unique += 1
            elif unique_data_alpha_ind_diff_corresponding_to_unique_event_id[i] == 1:
                ind_diff_1 += 1
                len_unique += 1

        if verbose == 1:
            print("% d'event sans diffusion = ", ind_diff_0 / len_unique)
            print("% d'event avec diffusion = ", ind_diff_1 / len_unique)

    ####################### Modification des ID de CPOP ###################################

    ################ data_alpha #########################################
    for ind_modif_id in range(0, len(data_alpha)):
        index_id_cell = np.where(real_id_cells == data_alpha[ind_modif_id]["ID_Cell"])
        data_alpha[ind_modif_id]["ID_Cell"] = perfect_id_cells[index_id_cell]

        index_cellule_emission = np.where(real_id_cells == data_alpha[ind_modif_id]["Cellule_D_Emission"])
        data_alpha[ind_modif_id]["Cellule_D_Emission"] = perfect_id_cells[index_cellule_emission]

    ################ data_edep_cell ######################################

    if test_file_not_empty != 0:
        elements_to_remove = []
        for ind_modif_id in range(0, len(data_edep_cell)):
            if ((data_edep_cell)[ind_modif_id]["ID_Cell"]) in deleted_id_txt:
                elements_to_remove.append(ind_modif_id)
        data_edep_cell = np.delete(data_edep_cell, elements_to_remove, 0)

    for ind_modif_id in range(0, len(data_edep_cell)):
        index_id_cell = np.where(real_id_cells == (data_edep_cell)[ind_modif_id]["ID_Cell"])
        (data_edep_cell)[ind_modif_id]["ID_Cell"] = perfect_id_cells[index_id_cell]

    ei = data_alpha["Ei"]  # Energy in keV
    ef = data_alpha["Ef"]

    n1 = number_of_lethal_events_for_alpha_traversals(dn1_de_continous_pre_calculated)

    n_tab = (n1(ei) - n1(ef))

    # print("############################# Calcul de dose #######################################")

    dosen_append_sur_une_simu_np = (((data_edep_cell)["fEdepn"]) * KEV_IN_J / masses_nuclei)
    dosec_append_sur_une_simu_np = (((data_edep_cell)["fEdepc"]) * KEV_IN_J / masses_cytoplasms)

    # test_dosen_append_sur_une_simu_np = np.sum(data_edep_cell[ind_dose]["fEdepn"], axis = 0)

    nb_particles_per_nucleus = np.bincount((data_alpha["ID_Cell"]).astype(int))
    ei_ef_unique_sur_une_simu  = np.bincount(data_alpha["ID_Cell"].astype(int), weights=ei - ef)

    while len(ei_ef_unique) < nb_cellules_reel:
        ei_ef_unique = np.append(ei_ef_unique, 0)

    while len(nb_particles_per_nucleus) < nb_cellules_reel:
        nb_particles_per_nucleus = np.append(nb_particles_per_nucleus, 0)

    #################################### Cross-fire dose au noyau ##########################################

    ############################################

    ind_non_cross_fire = data_alpha["ID_Cell"] == data_alpha["Cellule_D_Emission"]
    ind_cross_fire = data_alpha["ID_Cell"] != data_alpha["Cellule_D_Emission"]

    if indice_available_diffusion_info == 1:
        ind_non_cross_fire = ((data_alpha["ID_Cell"] == data_alpha["Cellule_D_Emission"]) &
                              (data_alpha["indice_if_diffusion"] == 0))
        ind_cross_fire = ((data_alpha["ID_Cell"] != data_alpha["Cellule_D_Emission"]) &
                          (data_alpha["indice_if_diffusion"] == 1))

    data_noyau_non_cross_fire = data_alpha[ind_non_cross_fire]

    dose_noyau_non_cross_fire = data_noyau_non_cross_fire["Ei"] - data_noyau_non_cross_fire["Ef"]

    dose_noyau_cross_fire = data_alpha["Ei"] - data_alpha["Ef"]

    dose_noyau_cross_fire = np.setdiff1d(dose_noyau_cross_fire, dose_noyau_non_cross_fire)

    sum_dose_noyau_crossfire = np.sum(dose_noyau_cross_fire)
    sum_dose_noyau_non_cross_fire = np.sum(dose_noyau_non_cross_fire)

    ################################################
    indice_zone1 = zone_cell == 1
    indice_zone2 = zone_cell == 2

    ei_ef_unique_non_cross_fire = np.bincount(((data_alpha["ID_Cell"])[ind_non_cross_fire]).astype(int),
                                              weights=ei[ind_non_cross_fire] - ef[ind_non_cross_fire])
    ei_ef_unique_cross_fire = np.bincount(((data_alpha["ID_Cell"])[ind_cross_fire]).astype(int),
                                          weights=ei[ind_cross_fire] - ef[ind_cross_fire])

    while len(ei_ef_unique_non_cross_fire) < nb_cellules_reel:
        ei_ef_unique_non_cross_fire = np.append(ei_ef_unique_non_cross_fire, 0)

    while len(ei_ef_unique_cross_fire) < nb_cellules_reel:
        ei_ef_unique_cross_fire = np.append(ei_ef_unique_cross_fire, 0)

    ei_ef_unique_non_cross_fire_zone1 = ei_ef_unique_non_cross_fire[indice_zone1]
    ei_ef_unique_non_cross_fire_zone2 = ei_ef_unique_non_cross_fire[indice_zone2]

    ei_ef_unique_cross_fire_zone1 = ei_ef_unique_cross_fire[indice_zone1]
    ei_ef_unique_cross_fire_zone2 = ei_ef_unique_cross_fire[indice_zone2]

    sum_dose_noyau_non_cross_fire_zone1 = np.sum(ei_ef_unique_non_cross_fire_zone1)
    sum_dose_noyau_non_cross_fire_zone2 = np.sum(ei_ef_unique_non_cross_fire_zone2)

    sum_dose_noyau_crossfire_zone1 = np.sum(ei_ef_unique_cross_fire_zone1)
    sum_dose_noyau_crossfire_zone2 = np.sum(ei_ef_unique_cross_fire_zone2)

    #################################### Nombre de cellules traversées par particule ##########################

    count_event_id = np.bincount((data_alpha["eventID"]).astype(int))

    for id_part in range(0, len(np.unique(data_alpha["eventID"]))):
        nb_nucl_traversees_par_la_particule = count_event_id[id_part]
        nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus.append(nb_nucl_traversees_par_la_particule)

    # print("############################# Calcul de survie #######################################")

    n_unique = np.bincount(data_alpha["ID_Cell"].astype(int), weights=n_tab)
    while len(n_unique) < nb_cellules_reel:
        n_unique = np.append(n_unique, 0)

    n_unique_tot_sur_une_simu = n_unique

    progress_bar['value'] += round(100 / (nb_complete_simulations - nb_files_with_errors), 2)
    progress_bar_label['text'] = update_progress_bar_label()
    window.update_idletasks()

    sum_dose_noyau_tot = sum_dose_noyau_crossfire + sum_dose_noyau_non_cross_fire
    analysis_dataframe['cross_fire_nucleus'] = sum_dose_noyau_crossfire / sum_dose_noyau_tot

    analysis_dataframe['ei_ef_sum'] = ei_ef_unique_sur_une_simu
    analysis_dataframe['nb_particles_per_nucleus'] = nb_particles_per_nucleus

    surviel_append_sur_une_simu = np.exp(-n_unique_tot_sur_une_simu)

    surviel_append_sur_une_simu[np.where(surviel_append_sur_une_simu == 0)] = 10 ** (-299)

    analysis_dataframe['cell_survival_local'] = surviel_append_sur_une_simu

    alpha_ref = 0.313  # HSG
    beta_ref  = 0.0615  # HSG

    dose_bio_append_sur_une_simu = \
        (np.sqrt(alpha_ref ** 2 - 4 * beta_ref * np.log(surviel_append_sur_une_simu)) - alpha_ref) / (2 * beta_ref)

    analysis_dataframe['biological_dose'] = dose_bio_append_sur_une_simu

    exp_surviel = np.exp(-np.asarray(surviel_append_sur_une_simu))
    tcp_une_simu = np.prod(exp_surviel)
    tcp_test_formula = np.prod(1 - surviel_append_sur_une_simu)
    analysis_dataframe['tcp'] = tcp_une_simu
    analysis_dataframe['tcp_test_formula'] = tcp_test_formula

    survieg_append_sur_une_simu = \
        np.exp(-n_unique_tot_sur_une_simu - BETAG[type_cell] * (dosen_append_sur_une_simu_np ** 2))

    analysis_dataframe['cell_survival_global'] = survieg_append_sur_une_simu

    analysis_dataframe['dose_nucleus'] = dosen_append_sur_une_simu_np
    analysis_dataframe['dose_cytoplasm'] = dosec_append_sur_une_simu_np
    analysis_dataframe['dose_cell'] = dosen_append_sur_une_simu_np + dosec_append_sur_une_simu_np

    spheroid_dose = data_edep_cell[0]["fEdep_sph"] * KEV_IN_J / masse_tum

    ################################################################################################################



# def todo_histo():
#     histo_nb_noy_par_p = 0 # mettre 1 affiche l'histogramme du nombre de noyaux traversés par les particules
#     histo_edep_noy_par_p = 0 # mettre 1 affiche l'histogramme de l'énergie moyenne déposée par particule dans un noyau
#                              # quand elle y rentre
#
#     if (histo_nb_noy_par_p==1) :
#         _, _, patches = plt.hist(nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus_np, bins=100,
#         edgecolor='white')
#
#         plt.xticks(fontsize=13.5)
#         plt.yticks(fontsize=13.5)
#         plt.xlabel('Nb of nuclei crossed by a particle',fontsize=15,fontname="Liberation Sans",fontweight='bold')
#         plt.ylabel('Occurence',fontsize=15,fontname="Liberation Sans",fontweight='bold')
#         #plt.title('Histogram of : Nb of nuclei crossed by a particle')
#         plt.grid(True)
#         plt.show()
#
#     if (histo_edep_noy_par_p==1) :
#         edep_in_nucleus_per_particle_np_histo=np.resize(edep_in_nucleus_per_particle_np,
#                   (1,len(edep_in_nucleus_per_particle_np)*len(edep_in_nucleus_per_particle_np[0])))
#         _, _, patches = plt.hist(edep_in_nucleus_per_particle_np_histo[0], bins=1000, edgecolor='black')
#
#         plt.xticks(fontsize=13.5)
#         plt.yticks(fontsize=13.5)
#         plt.xlabel('Edep by a particle in cell nucleus (keV)',fontsize=15,fontname="Liberation Sans",
#                   fontweight='bold')
#         plt.ylabel('Occurence',fontsize=15,fontname="Liberation Sans",fontweight='bold')
#         # plt.title('Histogram of : Edep by a particle in cell nucleus')
#         plt.grid(True)
#         plt.show()

def if_internalization_study():
    if labeling_percentage.winfo_exists():
        labeling_percentage.destroy()
    if labeling_percentage_entry.winfo_exists():
        labeling_percentage_entry.destroy()
    cell_compartment_label = tkinter.Label(window, text="Intra cellular distribution name :", fg='blue')
    cell_compartment_label.place(x=100, y=100)
    selected_distrib_name = tkinter.StringVar()
    cell_compartment_combobox = tkinter.ttk.Combobox(window, width=35, textvariable=selected_distrib_name)
    cell_compartment_combobox['values'] = ['Membrane', 'Cytoplasm', 'Homogeneous', 'Nucleus']
    cell_compartment_combobox.place(x=400, y=100)

def if_labeling_study() :
    global labeling_percentage,labeling_percentage_entry
    if cell_compartment_label.winfo_exists():
        cell_compartment_label.destroy()
    if cell_compartment_combobox.winfo_exists():
        cell_compartment_combobox.destroy()
    labeling_percentage = tkinter.Label(window, text="Labeling percentage : ", fg='blue')
    labeling_percentage.place(x=100, y=100)
    labeling_percentage_entry = tkinter.Entry(window, width=35)
    labeling_percentage_entry.place(x=400, y=100)

def id_deletion_of_root_outputs_with_errors():
    """
    Returns
    ======
        ids of simulations that returned a root output without errors when they were opened, in numpy format
        the number of errors encountered
    """
    indexe_of_root_output = 0
    indexes_root_files_without_errors = np.array([])

    start_time_test = time.perf_counter()

    while indexe_of_root_output < nb_complete_simulations:
        try: # Check if simulation 8 is valid for 10ppc (weird : there is diffusion)
            root_file_name =\
                f"Root/outputMultiCellulaire/{dossier_root}{nom_fichier_root}{indexe_of_root_output}_t0.root"
            with uproot.open(root_file_name) as root_file:
                _ = root_file['cell']['nameParticle'].array(library="np")
            indexes_root_files_without_errors = np.append(indexes_root_files_without_errors, indexe_of_root_output)
            indexe_of_root_output += 1
        except uproot.KeyInFileError:
            indexe_of_root_output += 1
    indexes_root_files_without_errors = np.sort(indexes_root_files_without_errors)
    nb_files_with_errors = nb_complete_simulations - len(indexes_root_files_without_errors)

    end_time_test = time.perf_counter()
    print("time while : ", end_time_test - start_time_test)
    return indexes_root_files_without_errors.astype(int), nb_files_with_errors

def graphic_window():
    global window, study_type_radiovalue, cell_compartment_combobox, geom_name_combobox, radionuclide_entry,\
        nb_simulations_entry, diffusion_combobox, number_particles_per_cell_combobox, cell_line_combobox,\
        cell_compartment_label, verbose_radiovalue

    window = tkinter.Tk()
    window.geometry("1000x700")

    study_type_label = tkinter.Label(window, text = "Type of study :", fg='red')
    study_type_label.place (x=100, y=50)

    labeling_percentage = tkinter.Label(window, text="Labeling percentage : ", fg='blue')
    labeling_percentage_entry = tkinter.Entry(window, width=35)

    study_type_radiovalue = tkinter.IntVar()
    study_type_radiovalue.set(0)
    study_type_radiovalue_0=tkinter.Radiobutton(window, text="Internalization",
                           variable=study_type_radiovalue,value=0, command=if_internalization_study)
    study_type_radiovalue_1=tkinter.Radiobutton(window, text="Labeling", variable=study_type_radiovalue,value=1,
                           command=if_labeling_study)
    study_type_radiovalue_0.place(x=390,y=50)
    study_type_radiovalue_1.place(x=590, y=50)

    cell_compartment_label = tkinter.Label(window, text="Intra cellular distribution name :", fg='blue')
    cell_compartment_label.place(x=100, y=100)
    selected_distrib_name = tkinter.StringVar()
    cell_compartment_combobox = tkinter.ttk.Combobox(window, width=35 , textvariable=selected_distrib_name)
    cell_compartment_combobox['values'] = ['Membrane', 'Cytoplasm', 'Homogeneous', 'Nucleus']
    cell_compartment_combobox.current(0)
    cell_compartment_combobox.place(x=400, y=100)

    geom_name_label = tkinter.Label(window, text = "Geometry name :", fg='blue')
    geom_name_label.place (x=100, y=150)
    geom_choice = tkinter.StringVar()
    geom_name_combobox = tkinter.ttk.Combobox(window, width=35 , textvariable=geom_choice)
    geom_name_combobox['values'] = ["30µmRadius Spheroid, 75 % cell packing", "50µmRadius Spheroid, 75 % cell packing",
                         "70µmRadius Spheroid, 75 % cell packing", "160µmRadius Spheroid, 75 % cell packing",
                         "95µmRadius Spheroid, 25 % cell packing", "95µmRadius Spheroid, 50 % cell packing",
                         "95µmRadius Spheroid, 75 % cell packing", "95µmRadius Spheroid, 75 % cell packing 2",
                         "100µmRadius Spheroid, 40 % cell packing"]
    geom_name_combobox.place(x=400, y=150)

    radionuclide_label = tkinter.Label(window, text = "Radionuclide used :", fg='blue')
    radionuclide_label.place (x=100, y=200)
    radionuclide_entry = tkinter.Entry(window, width=35)
    radionuclide_entry.insert(tkinter.END, "At211")
    radionuclide_entry.place(x=400, y=200)

    nb_simulations_label = tkinter.Label(window, text = "Number of simulations to analyse :", fg='blue')
    nb_simulations_label.place (x=100, y=250)
    nb_simulations_entry = tkinter.Entry(window, width=35)
    nb_simulations_entry.insert(tkinter.END, "20")
    nb_simulations_entry.place(x=400, y=250)

    diffusion_label = tkinter.Label(window, text = "Daughter diffusion :", fg='blue')
    diffusion_label.place (x=100, y=300)
    diffusion_choice = tkinter.StringVar()
    diffusion_combobox = tkinter.ttk.Combobox(window, width=35, textvariable=diffusion_choice)
    diffusion_combobox['values'] = ["Yes", "No"]
    diffusion_combobox.current(1)
    diffusion_combobox.place(x=400, y=300)

    number_particles_per_cell_label = tkinter.Label(window, text = "Number of alpha particles per cell :", fg='blue')
    number_particles_per_cell_label.place(x=100, y=350)
    number_particles_per_cell_choice = tkinter.StringVar()
    number_particles_per_cell_combobox = tkinter.ttk.Combobox(window, width=35,
                                                          textvariable=number_particles_per_cell_choice)
    number_particles_per_cell_combobox['values'] = ["1", "2", "3", "4", "5" ,"6", "7", "8", "9" ,"10", "42"]
    number_particles_per_cell_combobox.current(0)
    number_particles_per_cell_combobox.place(x=400, y=350)

    cell_line_label = tkinter.Label(window, text="Cell line :", fg='blue')
    cell_line_label.place(x=100, y=400)
    cell_line_choice = tkinter.StringVar()
    cell_line_combobox = tkinter.ttk.Combobox(window, width=35 , textvariable=cell_line_choice)
    cell_line_combobox['values'] = ['HSG', 'V79', 'CHO-K1']
    cell_line_combobox.current(0)
    cell_line_combobox.place(x=400, y=400)

    verbose_label = tkinter.Label(window, text = "Verbose :", fg='green')
    verbose_label.place (x=825, y=200)
    verbose_radiovalue=tkinter.IntVar()
    verbose_radiovalue.set(0)
    verbose_radiovalue_no=tkinter.Radiobutton(window, text="No", variable=verbose_radiovalue,value=0)
    verbose_radiovalue_yes=tkinter.Radiobutton(window, text="Yes", variable=verbose_radiovalue,value=1)
    verbose_radiovalue_yes.place(x=775,y=250)
    verbose_radiovalue_no.place(x=875, y=250)

    validate_button_1 = tkinter.Button(window, text = "Validate", command = add_new_buttons_to_graphic_window)
    validate_button_1.place(x=500, y=450)

    window.mainloop()

def create_folder_for_output_analysis_files():
    global dossier_root, index_of_first_root_output, nom_dossier_pour_excel_analyse
    dossier_root = study_type_folder_name + "/" + available_data_name_file[available_data_combobox.current()] + "/"
    index_of_first_root_output = 0 #Works only if the indexes of root files start at 0
    nom_dossier_pour_excel_analyse = f"{available_data_date[available_data_combobox.current()]}" \
                                     f"__{spheroid_compaction}CP_" \
                                     f"{r_sph}um_" \
                                     f"{rn_name}_diff{bool_diff[diffusion_combobox.current()]}_" \
                                     f"{nb_particles_per_cell[number_particles_per_cell_combobox.current()]}ppc_" \
                                     f"{cell_line_combobox.get()}"

    print("nom_dossier_pour_excel_analyse : ", nom_dossier_pour_excel_analyse)

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

def print_geometry_informations():
    print("nb_cellules_reel : ", nb_cellules_reel, '\n'
          "cell_packing = ", sum(masses_cells)/masse_tum, '\n'
          "masse_tum = ", masse_tum, '\n'
          "nb_cell_zone_1", nb_cell_zone_1, '\n'
          "nb_cell_zone_2", nb_cell_zone_2, '\n')
    return None

# def print_results():
#     return None

def main():
    global nb_cellules_reel, masses_cells, masse_tum, nb_cell_zone_1, nb_cell_zone_2
    verbose = verbose_radiovalue.get()
    create_folder_for_output_analysis_files()
    dn1_de_continous_pre_calculated = dn1_de_continous()
    ##################### Gestion des ID de CPOP ##################################################
    txt_id_deleted_cells = f"Cpop_Deleted_Cells_ID_Txt/IDCell_{nom_config}.txt"
    real_id_cells, test_file_not_empty, deleted_id_txt =\
        geometry_informations.cpop_real_cell_id_determination(txt_id_deleted_cells, nb_cellules_xml)
    nb_cellules_reel = len(real_id_cells)
    perfect_id_cells = np.arange(0,nb_cellules_reel)
    ###################### Lecture Geométrie #####################################
    ###### Masses #######
    txt_cells_masses=f"Cpop_Masse_Txt/MassesCell_{nom_config}.txt"
    masses_cytoplasms, masses_nuclei, masses_cells = geometry_informations.masses_cells_reading(txt_cells_masses)
    r_tum = float(r_sph) * 10**(-6) #in meters
    masse_tum=((4/3)*np.pi*r_tum**3)*1000 #in kg
    ###### Positions ######
    positions_x, positions_y, positions_z = geometry_informations.positions_cells_reading(xml_geom, real_id_cells)

    zone_cell, nb_cell_zone_1, nb_cell_zone_2 = determine_cells_in_2_spheroid_zones(positions_x,
                                                    positions_y, positions_z,
                                                    radius_zone_1 = 50, radius_zone_2 = 95,
                                                    nb_cells = nb_cellules_reel)
    if verbose == 1:
        print_geometry_informations()

    ######################## Root ##################################################################
    nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus, surviel_append_sur_une_simu,\
        surviel_append_sur_toutes_simus, survieg_append_sur_une_simu, survieg_append_sur_toutes_simus,\
        dosen_append_sur_une_simu, dosen_append_sur_toutes_simus, dosec_append_sur_une_simu,\
        dosec_append_sur_toutes_simus, dosem_append_sur_une_simu,\
        dosen_c_append_sur_toutes_simus, ratio_crossfire_noyau_sur_toutes_simus,\
        ratio_crossfire_noyau_sur_toutes_simus_zone1, ratio_crossfire_noyau_sur_toutes_simus_zone2,\
        ei_ef_sum_sur_toutes_simus, nombre_particules_par_noyau_sur_toutes_simus, tcp_append_sur_toutes_simus,\
        tcp_test_formula_append_sur_toutes_simus, simulation_id, dose_bio_append_sur_toutes_simus, test_spectre_diff =\
        ([] for _ in range(21))

    indexes_root_files_without_errors_np, nb_files_with_errors = \
        id_deletion_of_root_outputs_with_errors()

    ##################################################

    for i in indexes_root_files_without_errors_np:
        sum_dose_noyau_crossfire = sum_dose_noyau_non_cross_fire = sum_dose_noyau_non_cross_fire_zone1 = \
        sum_dose_noyau_non_cross_fire_zone2 = sum_dose_noyau_crossfire_zone1 = sum_dose_noyau_crossfire_zone2 = 0

        dosen_append_sur_une_simu_np = np.zeros(nb_cellules_reel)
        dosec_append_sur_une_simu_np = np.zeros(nb_cellules_reel)
        dosem_append_sur_une_simu_np = np.zeros(nb_cellules_reel)

        n_unique_tot_sur_une_simu = np.zeros(nb_cellules_reel)
        ei_ef_unique_sur_une_simu = np.zeros(nb_cellules_reel)
        nb_particles_per_nucleus = np.zeros(nb_cellules_reel)

        simulation_id.append(i) #Works only if no multithreading.

        root_file_name = f"Root/outputMultiCellulaire/{dossier_root}{nom_fichier_root}{i}_t0.root"

        if verbose == 1:
            print("Root file name : ", root_file_name)

        root_data_opened, ind_alphaplusplus, ind_alphaplus, ind_helium, data_alpha, ind_end_of_run, data_edep_cell = \
            ([] for _ in range(7))
        indice_available_diffusion_info = indice_available_edep_sph_info = 0

        f = uproot.open(root_file_name)
        name_particle_root = f['cell']['nameParticle'].array(library="np")

        print("len(name_particle_root)", len(name_particle_root))
        ei_root = f['cell']['Ei'].array(library="np")
        print("ei_root", len(ei_root))
        ef_root = f['cell']['Ef'].array(library="np")
        id_cell_root = f['cell']['ID_Cell'].array(library="np")
        emission_cell_root = f['cell']['Cellule_D_Emission'].array(library="np")
        event_id_root = f['cell']['eventID'].array(library="np")
        energy_deposited_nucleus_root = f['cell']['fEdepn'].array(library="np")
        energy_deposited_cytoplasm_root = f['cell']['fEdepc'].array(library="np")
        try:
            diffusion_index_root = f['cell']['indice_if_diffusion'].array(library="np")
            indice_available_diffusion_info = 1
        except:
            print("indice_if_diffusion not available on these data")
            indice_available_diffusion_info = 0
        try:
            energy_deposited_spheroid_root = f['cell']['fEdep_sph'].array(library="np")
            indice_available_edep_sph_info = 1
        except:
            print("fEdep_sph not available on these data")
            indice_available_edep_sph_info = 0

        if indice_available_diffusion_info == 0:
            root_data_opened.append(np.core.records.fromarrays([name_particle_root, ei_root, ef_root,
                                                                id_cell_root, emission_cell_root,
                                                                event_id_root, energy_deposited_nucleus_root,
                                                                energy_deposited_cytoplasm_root],
                           names='nameParticle, Ei, Ef, ID_Cell, Cellule_D_Emission, eventID, fEdepn, fEdepc'))
        elif indice_available_diffusion_info == 1 and indice_available_edep_sph_info == 0:
            root_data_opened.append(np.core.records.fromarrays([name_particle_root, ei_root, ef_root,
                                                                id_cell_root, emission_cell_root,
                                                                event_id_root, energy_deposited_nucleus_root,
                                                                energy_deposited_cytoplasm_root,
                                                                diffusion_index_root],
                        names='nameParticle, Ei, Ef, ID_Cell, Cellule_D_Emission, eventID, fEdepn, fEdepc,'
                              ' indice_if_diffusion'))
        elif indice_available_diffusion_info == 1 and indice_available_edep_sph_info == 1:
            root_data_opened.append(np.core.records.fromarrays([name_particle_root, ei_root, ef_root,
                                                                id_cell_root, emission_cell_root,
                                                                event_id_root, energy_deposited_nucleus_root,
                                                                energy_deposited_cytoplasm_root,
                                                                diffusion_index_root,
                                                                energy_deposited_spheroid_root],
                        names='nameParticle, Ei, Ef, ID_Cell, Cellule_D_Emission, eventID, fEdepn, fEdepc,'
                              ' indice_if_diffusion, fEdep_sph'))
            print("len(root_data_opened[0])", len(root_data_opened[0]))
            print("root_data_opened[0][:100]", root_data_opened[0][:100])

            ind_division_simus = 0

        ind_alphaplusplus.append((root_data_opened[ind_division_simus])["nameParticle"]=='alpha') #Check format (remove ind_division_simus)
        ind_alphaplus.append((root_data_opened[ind_division_simus])["nameParticle"]=='alpha+')
        ind_helium.append((root_data_opened[ind_division_simus])["nameParticle"]=='helium')

        data_alpha.append((np.concatenate(((root_data_opened[ind_division_simus])[ind_alphaplusplus[ind_division_simus]],
                        (root_data_opened[ind_division_simus])[ind_alphaplus[ind_division_simus]],
                        (root_data_opened[ind_division_simus])[ind_helium[ind_division_simus]]))))

        (data_alpha[ind_division_simus])["eventID"] +=\
            (nb_cellules_reel / nb_group_of_cells_considered)*ind_division_simus

        ind_end_of_run.append((root_data_opened[ind_division_simus])["nameParticle"] == 'EndOfRun')

        data_edep_cell.append((root_data_opened[ind_division_simus])[ind_end_of_run[ind_division_simus]])

        data_alpha = np.concatenate([row for row in data_alpha])

        ########################## Vérification diffusion aux bonnes énergies ###############################

        if indice_available_diffusion_info == 1:

            unique_data_alpha_event_id = np.unique(data_alpha['eventID'], return_index=True)

            ind_diff_0 = ind_diff_1 = len_unique = 0

            indices_ab = unique_data_alpha_event_id[1]

            unique_data_alpha_ind_diff_corresponding_to_unique_event_id =\
                np.take(data_alpha['indice_if_diffusion'], indices_ab)

            for i in range(0,len(unique_data_alpha_ind_diff_corresponding_to_unique_event_id)):
                if unique_data_alpha_ind_diff_corresponding_to_unique_event_id[i] == 0:
                    ind_diff_0+=1
                    len_unique+=1
                elif unique_data_alpha_ind_diff_corresponding_to_unique_event_id[i] == 1:
                    ind_diff_1 += 1
                    len_unique += 1

            if verbose == 1:
                print("% d'event sans diffusion = ", ind_diff_0/len_unique)
                print("% d'event avec diffusion = ", ind_diff_1/len_unique)

        ####################### Modification des ID de CPOP ###################################

        ################ data_alpha #########################################
        for ind_modif_id in range(0,len(data_alpha)):
            index_id_cell = np.where(real_id_cells == data_alpha[ind_modif_id]["ID_Cell"])
            data_alpha[ind_modif_id]["ID_Cell"] = perfect_id_cells[index_id_cell]

            index_cellule_emission = np.where(real_id_cells == data_alpha[ind_modif_id]["Cellule_D_Emission"])
            data_alpha[ind_modif_id]["Cellule_D_Emission"] = perfect_id_cells[index_cellule_emission]


        ################ data_edep_cell ######################################

        if test_file_not_empty != 0:

            for ind_dose in range(0, nb_group_of_cells_considered):
                elements_to_remove = []
                for ind_modif_id in range(0, len(data_edep_cell[ind_dose])):
                    if ((data_edep_cell[ind_dose])[ind_modif_id]["ID_Cell"]) in deleted_id_txt:
                        elements_to_remove.append(ind_modif_id)
                data_edep_cell[ind_dose] = np.delete(data_edep_cell[ind_dose], elements_to_remove, 0)

        for ind_dose in range(0, nb_group_of_cells_considered):
            for ind_modif_id in range(0,len(data_edep_cell[ind_dose])):
                index_id_cell = np.where(real_id_cells == (data_edep_cell[ind_dose])[ind_modif_id]["ID_Cell"])
                (data_edep_cell[ind_dose])[ind_modif_id]["ID_Cell"] = perfect_id_cells[index_id_cell]

        ei = data_alpha["Ei"] # Energy in keV
        ef = data_alpha["Ef"]

        n1 = number_of_lethal_events_for_alpha_traversals(dn1_de_continous_pre_calculated)

        n_tab = (n1(ei) - n1(ef))

        # print("############################# Calcul de dose #######################################")


        for ind_dose in range(0, nb_group_of_cells_considered):
            dosen_append_sur_une_simu_np += (((data_edep_cell[ind_dose])["fEdepn"]) * KEV_IN_J / masses_nuclei)
            dosec_append_sur_une_simu_np += (((data_edep_cell[ind_dose])["fEdepc"]) * KEV_IN_J / masses_cytoplasms)

        # test_dosen_append_sur_une_simu_np = np.sum(data_edep_cell[ind_dose]["fEdepn"], axis = 0)

        count_cell_id = np.bincount((data_alpha["ID_Cell"]).astype(int))
        ei_ef_unique = np.bincount(data_alpha["ID_Cell"].astype(int), weights= ei - ef)

        while len(ei_ef_unique) < nb_cellules_reel:
            ei_ef_unique = np.append(ei_ef_unique, 0)

        while len(count_cell_id) < nb_cellules_reel:
            count_cell_id = np.append(count_cell_id, 0)

        ei_ef_unique_sur_une_simu += ei_ef_unique
        nb_particles_per_nucleus += count_cell_id


        #################################### Cross-fire dose au noyau ##########################################


        ############################################

        ind_non_cross_fire = data_alpha["ID_Cell"] == data_alpha["Cellule_D_Emission"]
        ind_cross_fire = data_alpha["ID_Cell"] != data_alpha["Cellule_D_Emission"]

        if indice_available_diffusion_info == 1:
            ind_non_cross_fire = ((data_alpha["ID_Cell"] == data_alpha["Cellule_D_Emission"]) &
                                (data_alpha["indice_if_diffusion"]==0))
            ind_cross_fire = ((data_alpha["ID_Cell"] != data_alpha["Cellule_D_Emission"]) &
                             (data_alpha["indice_if_diffusion"]==1))

        data_noyau_non_cross_fire = data_alpha[ind_non_cross_fire]

        dose_noyau_non_cross_fire=data_noyau_non_cross_fire["Ei"]-data_noyau_non_cross_fire["Ef"]

        dose_noyau_cross_fire=data_alpha["Ei"]-data_alpha["Ef"]

        dose_noyau_cross_fire=np.setdiff1d(dose_noyau_cross_fire,dose_noyau_non_cross_fire)

        sum_dose_noyau_crossfire+=np.sum(dose_noyau_cross_fire)
        sum_dose_noyau_non_cross_fire+=np.sum(dose_noyau_non_cross_fire)

        ################################################
        indice_zone1 = zone_cell == 1
        indice_zone2 = zone_cell == 2

        ei_ef_unique_non_cross_fire = np.bincount(((data_alpha["ID_Cell"])[ind_non_cross_fire]).astype(int),
                                                weights=ei[ind_non_cross_fire] - ef[ind_non_cross_fire])
        ei_ef_unique_cross_fire = np.bincount(((data_alpha["ID_Cell"])[ind_cross_fire]).astype(int),
                                                weights=ei[ind_cross_fire] - ef[ind_cross_fire])

        while len(ei_ef_unique_non_cross_fire) < nb_cellules_reel:
            ei_ef_unique_non_cross_fire = np.append(ei_ef_unique_non_cross_fire, 0)

        while len(ei_ef_unique_cross_fire) < nb_cellules_reel:
            ei_ef_unique_cross_fire = np.append(ei_ef_unique_cross_fire, 0)

        ei_ef_unique_non_cross_fire_zone1 = ei_ef_unique_non_cross_fire[indice_zone1]
        ei_ef_unique_non_cross_fire_zone2 = ei_ef_unique_non_cross_fire[indice_zone2]

        ei_ef_unique_cross_fire_zone1 = ei_ef_unique_cross_fire[indice_zone1]
        ei_ef_unique_cross_fire_zone2 = ei_ef_unique_cross_fire[indice_zone2]

        sum_dose_noyau_non_cross_fire_zone1 += np.sum(ei_ef_unique_non_cross_fire_zone1)
        sum_dose_noyau_non_cross_fire_zone2 += np.sum(ei_ef_unique_non_cross_fire_zone2)

        sum_dose_noyau_crossfire_zone1 += np.sum(ei_ef_unique_cross_fire_zone1)
        sum_dose_noyau_crossfire_zone2 += np.sum(ei_ef_unique_cross_fire_zone2)


        #################################### Nombre de cellules traversées par particule ##########################

        count_event_id = np.bincount((data_alpha["eventID"]).astype(int))

        for id_part in range(0, len(np.unique(data_alpha["eventID"]))):
            nb_nucl_traversees_par_la_particule=count_event_id[id_part]
            nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus.append(nb_nucl_traversees_par_la_particule)


        # print("############################# Calcul de survie #######################################")

        n_unique=np.bincount(data_alpha["ID_Cell"].astype(int), weights=n_tab)
        while len(n_unique) < nb_cellules_reel:
            n_unique=np.append(n_unique,0)

        n_unique_tot_sur_une_simu+=n_unique

        progress_bar['value'] += round(100 / (nb_complete_simulations - nb_files_with_errors), 2)
        progress_bar_label['text'] = update_progress_bar_label()
        window.update_idletasks()

        sum_dose_noyau_tot=sum_dose_noyau_crossfire+sum_dose_noyau_non_cross_fire
        ratio_crossfire_noyau_sur_une_simu=sum_dose_noyau_crossfire/sum_dose_noyau_tot
        ratio_crossfire_noyau_sur_toutes_simus.append(ratio_crossfire_noyau_sur_une_simu)

        ratio_crossfire_noyau_sur_une_simu_zone1=0
        ratio_crossfire_noyau_sur_toutes_simus_zone1.append(ratio_crossfire_noyau_sur_une_simu_zone1)

        ratio_crossfire_noyau_sur_une_simu_zone2=0
        ratio_crossfire_noyau_sur_toutes_simus_zone2.append(ratio_crossfire_noyau_sur_une_simu_zone2)

        ei_ef_sum_sur_toutes_simus.append(ei_ef_unique_sur_une_simu)
        nombre_particules_par_noyau_sur_toutes_simus.append(nb_particles_per_nucleus)

        surviel_append_sur_une_simu=np.exp(-n_unique_tot_sur_une_simu)

        surviel_append_sur_une_simu[np.where(surviel_append_sur_une_simu == 0)] = 10**(-299)

        surviel_append_sur_toutes_simus.append(surviel_append_sur_une_simu)

        alpha_ref = 0.313  # HSG
        beta_ref = 0.0615  # HSG

        dose_bio_append_sur_une_simu =\
            (np.sqrt(alpha_ref ** 2 - 4 * beta_ref * np.log(surviel_append_sur_une_simu)) - alpha_ref) / (2 * beta_ref)
        dose_bio_append_sur_toutes_simus.append(dose_bio_append_sur_une_simu)

        exp_surviel = np.exp(-np.asarray(surviel_append_sur_une_simu))
        tcp_une_simu = np.prod(exp_surviel)
        tcp_test_formula = np.prod(1 - surviel_append_sur_une_simu)
        tcp_append_sur_toutes_simus.append(tcp_une_simu)
        tcp_test_formula_append_sur_toutes_simus.append(tcp_test_formula)

        survieg_append_sur_une_simu = \
            np.exp(-n_unique_tot_sur_une_simu - BETAG[type_cell] * (dosen_append_sur_une_simu_np ** 2))
        survieg_append_sur_toutes_simus.append(survieg_append_sur_une_simu)


        dosen_append_sur_toutes_simus.append(dosen_append_sur_une_simu_np)
        dosec_append_sur_toutes_simus.append(dosec_append_sur_une_simu_np)
        dosen_c_append_sur_toutes_simus.append(dosen_append_sur_une_simu_np+dosec_append_sur_une_simu_np)

        spheroid_dose = data_edep_cell[0]["fEdep_sph"]* KEV_IN_J / masse_tum

        ################################################################################################################

    progress_bar['value'] = math.floor(progress_bar['value'])

    surviel_append_sur_toutes_simus_np=np.asarray(surviel_append_sur_toutes_simus)
    survieg_append_sur_toutes_simus_np=np.asarray(survieg_append_sur_toutes_simus)


    ei_ef_sum_sur_toutes_simus_np = np.asarray(ei_ef_sum_sur_toutes_simus)
    nombre_particules_par_noyau_sur_toutes_simus_np = np.asarray(nombre_particules_par_noyau_sur_toutes_simus)
    nombre_particules_par_noyau_sur_toutes_simus_sans_zero_np = np.array([])
    ei_ef_sum_sur_toutes_simus_sans_zero_np = np.array([])

    #Suppression des cells qui n'ont pas été traversées par des alpha,
    # pour le calcul d'Edep moy quand les particules touchent des noyaux ###########

    nb_cells_suppr = 0

    for i in range(0,(nb_complete_simulations - nb_files_with_errors)):
        ind_nucleus_with_no_particles = [index for index, value in enumerate(
            nombre_particules_par_noyau_sur_toutes_simus_np[i]) if value == 0]

        nombre_particules_par_noyau_sur_toutes_simus_sans_zero_np =\
            np.append(nombre_particules_par_noyau_sur_toutes_simus_sans_zero_np,
                      np.delete(nombre_particules_par_noyau_sur_toutes_simus_np[i], ind_nucleus_with_no_particles))
        ei_ef_sum_sur_toutes_simus_sans_zero_np =\
            np.append(ei_ef_sum_sur_toutes_simus_sans_zero_np,
                      np.delete(ei_ef_sum_sur_toutes_simus_np[i], ind_nucleus_with_no_particles))

        nb_cells_suppr += len(ind_nucleus_with_no_particles)


    print()

    edep_in_nucleus_per_particle_np = \
        ei_ef_sum_sur_toutes_simus_sans_zero_np/nombre_particules_par_noyau_sur_toutes_simus_sans_zero_np
    mean_edep_dans_noy_par_particule = np.mean(edep_in_nucleus_per_particle_np)

    dosen_append_sur_toutes_simus_np=np.asarray(dosen_append_sur_toutes_simus)
    dosec_append_sur_toutes_simus_np=np.asarray(dosec_append_sur_toutes_simus)
    dosen_c_append_sur_toutes_simus_np=np.asarray(dosen_c_append_sur_toutes_simus)
    ratio_crossfire_noyau_sur_toutes_simus_np=np.asarray(ratio_crossfire_noyau_sur_toutes_simus)
    ratio_crossfire_noyau_sur_toutes_simus_zone1_np=np.asarray(ratio_crossfire_noyau_sur_toutes_simus_zone1)
    ratio_crossfire_noyau_sur_toutes_simus_zone2_np=np.asarray(ratio_crossfire_noyau_sur_toutes_simus_zone2)
    tcp_append_sur_toutes_simus_np = np.asarray(tcp_append_sur_toutes_simus)
    tcp_test_formula_append_sur_toutes_simus_np = np.asarray(tcp_test_formula_append_sur_toutes_simus)
    dose_bio_append_sur_toutes_simus_np = np.asarray(dose_bio_append_sur_toutes_simus)

    surviel=np.mean(surviel_append_sur_toutes_simus_np, axis=0)
    survieg=np.mean(survieg_append_sur_toutes_simus_np, axis=0)
    dosen_tot=np.mean(dosen_append_sur_toutes_simus_np, axis=0)
    dosec_tot=np.mean(dosec_append_sur_toutes_simus_np, axis=0)
    dosen_c_tot=np.mean(dosen_c_append_sur_toutes_simus_np, axis=0)
    tcp = np.mean(tcp_append_sur_toutes_simus_np, axis=0)
    tcp_formula2 = np.mean(tcp_test_formula_append_sur_toutes_simus_np,axis=0)
    dose_bio = np.mean(dose_bio_append_sur_toutes_simus_np, axis=0)
    mean_cross_fire = np.full(nb_cellules_reel ,np.mean(ratio_crossfire_noyau_sur_toutes_simus_np)*100)

    dose_tot=dosen_tot+dosec_tot

    sum_dose_noyau = np.mean(np.sum(dosen_append_sur_toutes_simus_np,axis=1))
    sum_dose_tot = np.mean(np.sum(dosen_c_append_sur_toutes_simus_np,axis=1))

    incertdmoy_n_tot=np.std(dosen_append_sur_toutes_simus_np,axis=0)
    incertdmoy_c_tot=np.std(dosec_append_sur_toutes_simus_np,axis=0)
    incertdmoy_n_c_tot=np.std(dosen_c_append_sur_toutes_simus_np,axis=0)
    incert_crossfire=np.full(nb_cellules_reel ,np.std(ratio_crossfire_noyau_sur_toutes_simus_np)*100)
    incert_mean_cell_survival=np.full(nb_cellules_reel ,np.std(np.mean(surviel_append_sur_toutes_simus_np,axis=1)))
    incert_dose_bio=np.std(np.mean(dose_bio_append_sur_toutes_simus_np, axis=1))

    alpha=0.5

    sf=np.sum(np.exp(-alpha*dosen_c_tot))/nb_cellules_reel

    eud=-np.log(sf)/alpha

    error_nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus = \
       np.std(nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus)\
            /np.sqrt(len(nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus))

    print(f"TCP = {tcp} +- {2*np.std(tcp_append_sur_toutes_simus_np)/np.sqrt(nb_complete_simulations)}", '\n'
          f"TCP test formula = {tcp_formula2} +- "
                      f"{2*np.std(tcp_test_formula_append_sur_toutes_simus_np)/np.sqrt(nb_complete_simulations)}", '\n'
          f"EUD = {eud}", '\n'
          f"Min dose totale absorbée = {np.mean(np.min(dosen_c_append_sur_toutes_simus_np,axis=1))} +- "
                      f"{2 * np.std(np.min(dosen_c_append_sur_toutes_simus_np, axis=1))}", '\n'
          f"Max dose totale absorbée = {np.mean(np.max(dosen_c_append_sur_toutes_simus_np,axis=1)):.2f} +- "
                      f"{(2 * np.std(np.max(dosen_c_append_sur_toutes_simus_np, axis=1))):.2f}", '\n'
          # Les 2 méthodes de calcul de la moyenne sont les mêmes, mais la seconde calcule correctement l'incertitude
          f"Moyenne des doses absorbées aux noyaux = {np.mean(np.mean(dosen_append_sur_toutes_simus_np,axis=1))} +- "
                      f"{2*np.std(np.mean(dosen_append_sur_toutes_simus_np,axis=1))}", '\n'
          f"Moyenne des doses absorbées aux cytoplasmes = {np.mean(np.mean(dosec_append_sur_toutes_simus_np,axis=1))}"
                      f" +- {2*np.std(np.mean(dosec_append_sur_toutes_simus_np,axis=1))}", '\n'
          f"Moyenne des doses absorbées aux cellules = {np.mean(np.mean(dosen_c_append_sur_toutes_simus_np,axis=1))}"
                      f" +- {2*np.std(np.mean(dosen_c_append_sur_toutes_simus_np,axis=1))}", '\n'
          f"Sum dose cellules = {np.sum(dosen_c_tot):.2f}", '\n'
          f"Energie moyenne déposée par une particule quand elle touche un noyau = {mean_edep_dans_noy_par_particule}"
                      f" +- {np.std(edep_in_nucleus_per_particle_np)/np.sqrt(len(edep_in_nucleus_per_particle_np))}",
          '\n'
          f"Nombre moyen de noyaux traversés par une particule en moyenne ="
                      f"{np.mean(nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus)} +-"
                      f"{error_nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus}", '\n'
          f"Cross-fire dose au noyau en moyenne = {np.mean(ratio_crossfire_noyau_sur_toutes_simus_np)*100} % +- "
                      f"{2*np.std(ratio_crossfire_noyau_sur_toutes_simus_np)*100} %", '\n'
          f"Cross-fire dose au noyau en moyenne = {np.mean(ratio_crossfire_noyau_sur_toutes_simus_zone1_np)*100} % +- "
                      f"{2*np.std(ratio_crossfire_noyau_sur_toutes_simus_zone1_np)*100} %", '\n'
          f"Cross-fire dose au noyau en moyenne = {np.mean(ratio_crossfire_noyau_sur_toutes_simus_zone2_np)*100} % +- "
                      f"{2*np.std(ratio_crossfire_noyau_sur_toutes_simus_zone2_np)*100} %", '\n'
          f"Moyenne des survies cellulaires = {np.mean(np.mean(surviel_append_sur_toutes_simus_np,axis=1))} +- "
                      f"{2*np.std(np.mean(surviel_append_sur_toutes_simus_np,axis=1))}", '\n'
          f"Moyenne des doses biologiques = {np.mean(np.mean(dose_bio_append_sur_toutes_simus_np, axis=1))} +- "
                      f"{2 * np.std(np.mean(dose_bio_append_sur_toutes_simus_np, axis=1))}", '\n'
          f"Max des doses biologiques = {np.mean(np.max(dose_bio_append_sur_toutes_simus_np, axis=1))} +- "
                      f"{2 * np.std(np.max(dose_bio_append_sur_toutes_simus_np, axis=1))}", '\n'
          )

    if indice_available_edep_sph_info:
        print("Dose sphéroïde")
        print(spheroid_dose[0], "Gy")

    print("Nombre de simus fonctionnelles = ")
    print(nb_complete_simulations - nb_files_with_errors)

    print()

    id_cell_arr=np.arange(nb_cellules_reel)


    fname = "AnalysisResults/" + study_type_folder_name + "/" + nom_dossier_pour_excel_analyse + "/" + "Emission" +\
            cell_compartment + ".xlsx"
    name_columns = ['ID_Cell','Zone_Cell' ,'Survie locale', 'Incert Survie Locale','Survie globale',
                    'Dmoy au noyau, par simu (Gy)','Incert Dmoy au noyau par simu(Gy)','Dmoy cyto, par simu (Gy)',
                    'incert Dmoy cyto par simu (Gy)','Dmoy au noyau+cyto, par simu (Gy)',
                    'incert Dmoy noy+cyto par simu (Gy)', 'Ratio Crossfire %','Incert Ratio Crossfire %',
                    'Ratio Crossfire zone 1 %', 'Ratio Crossfire zone 2 %','Nb noyau par particule',
                    'Sum dose sur tous noyaux', 'Sum dose sur toutes cellules', 'Edep moy par particule dans noyau',
                    'Incert Edep moy par particule dans noyau', 'Nb noyaux par particule',
                    'Incert nb noyaux par particule' ,'TCP', 'Incert_TCP','TCP formula 2','Incert TCP formula 2',
                    'EUD','Dose bio','IncertDoseBio', "DoseSpheroid"]
    results_to_write = np.array(np.concatenate([[name_columns],
                            np.transpose([id_cell_arr,zone_cell,surviel, incert_mean_cell_survival, survieg,
                                          dosen_tot,np.full(nb_cellules_reel,incertdmoy_n_tot),dosec_tot,
                                          np.full(nb_cellules_reel,incertdmoy_c_tot),dosen_c_tot,
                                          np.full(nb_cellules_reel,incertdmoy_n_c_tot), mean_cross_fire,
                                          incert_crossfire,np.full(nb_cellules_reel,0), np.full(nb_cellules_reel,0),
                                          np.full(nb_cellules_reel, np.full(nb_cellules_reel,
                                               np.mean(nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus))),
                                          np.full(nb_cellules_reel,sum_dose_noyau),
                                          np.full(nb_cellules_reel,sum_dose_tot),
                                          np.full(nb_cellules_reel,mean_edep_dans_noy_par_particule),
                                          np.full(nb_cellules_reel, np.std(edep_in_nucleus_per_particle_np)
                                               /np.sqrt(len(edep_in_nucleus_per_particle_np))),
                                          np.full(nb_cellules_reel,
                                               np.mean(nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus)),
                                          np.full(nb_cellules_reel,
                                               np.std(nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus)
                                               /np.sqrt(len(nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus))),
                                          np.full(nb_cellules_reel,tcp),
                                          np.full(nb_cellules_reel, np.std(tcp_append_sur_toutes_simus_np)),
                                          np.full(nb_cellules_reel,tcp_formula2),
                                          np.full(nb_cellules_reel,np.std(tcp_test_formula_append_sur_toutes_simus_np)),
                                          np.full(nb_cellules_reel,eud),
                                          dose_bio, np.full(nb_cellules_reel,incert_dose_bio),
                                          np.full(nb_cellules_reel,spheroid_dose[0])])]))
    wb = Workbook()
    wb.new_sheet("AnalysisResults", data=results_to_write)
    wb.save(fname)


    fname = "AnalysisResults/" + study_type_folder_name + "/" + nom_dossier_pour_excel_analyse +\
            "/" + "Emission" + cell_compartment + "SimID" + ".xlsx"
    name_columns = ['SimulationID', 'TCP', 'Dose moyenne noyau', 'Incert dose moyenne noyau', 'Dose moyenne cellule',
                    'Incert dose moyenne cellule', 'Dose bio moyenne', 'Incert dose bio moyenne']
    results_to_write =\
        np.array(np.concatenate([[name_columns],
             np.transpose([simulation_id, tcp_append_sur_toutes_simus_np,
                         np.mean(dosen_append_sur_toutes_simus_np,axis=1),
                         np.std(dosen_append_sur_toutes_simus_np,axis=1)/np.sqrt(nb_cellules_reel),
                         np.mean(dosen_c_append_sur_toutes_simus_np,axis=1),
                         np.std(dosen_c_append_sur_toutes_simus_np,axis=1)/np.sqrt(nb_cellules_reel),
                         np.mean(dose_bio_append_sur_toutes_simus_np,axis=1),
                         np.std(dose_bio_append_sur_toutes_simus_np,axis=1)/np.sqrt(nb_cellules_reel)])]))
    wb = Workbook()
    wb.new_sheet("AnalysisResults_SimID", data=results_to_write)
    wb.save(fname)

    # progress_bar.stop()
    progress_bar_label['text'] = update_progress_bar_label()

    end_time = time.perf_counter()

    print(" Temps total =  ", (end_time - start_time)//60, "minutes",
          math.floor((end_time - start_time)%60), "secondes")

def update_progress_bar_label():
    return f"Current progress: {progress_bar['value']} %"

def add_new_buttons_to_graphic_window():
    global r_sph, nom_config, spheroid_compaction, xml_geom, nb_cellules_xml, cell_compartment,\
    nb_complete_simulations, nb_group_of_cells_considered, simulation_name,\
    study_type_folder_name, bool_diff, rn_name, nb_particles_per_cell, type_cell, available_data_date,\
    available_data_name_file, available_data_combobox, nom_fichier_root, progress_bar, progress_bar_label

    geom_list = ["Elg030um75CP", "Elg050um75CP", "Elg070um75CP", "Elg160um75CP", "Elg095um25CP",
                 "Elg095um50CP", "Elg095um75CP", "Elg095um75CP_2", "Elg100um40CP"]
    r_sph = geom_list[geom_name_combobox.current()][3:6]

    nom_config = (geom_list[geom_name_combobox.current()])  # Les fichiers contenant les masses de toutes les cellules,
                                                 # et ceux des ID de cellules supprimés de CPOP à G4,
                                                 # sont appelés MassesCell_nom_config.txt, et IDCell_nom_config.txt
    spheroid_compaction = geom_list[geom_name_combobox.current()][8:10]

    xml_geom = "Cpop_Geom_XML/" + nom_config + ".cfg" + ".xml"

    nb_cellules_xml = geometry_informations.count_number_of_cells_in_xml_file(xml_geom)
    # Nombre de cellules contenues dans le fichier .xml de géométrie créé par CPOP

    cell_compartment = (cell_compartment_combobox.get())

    study_type = study_type_radiovalue.get()  # 0 for internalization study, 1 for labeling study

    nb_complete_simulations = int(nb_simulations_entry.get())
    nb_group_of_cells_considered = 1  #Usually 1 but can be changed for very long simulations.
                                      # E.g. : if = 2, half the jobs were sent in 50% of the cells
                                      #and the other jobs in the other 50% of cells.

    if study_type == 0:
        simulation_name = cell_compartment
        study_type_folder_name = "Internalization"
    elif study_type == 1:
        labeling_percentage_get = labeling_percentage_entry.get()
        labeling_percentage_name = str(labeling_percentage_get) + "_Percent"
        simulation_name = labeling_percentage_name
        study_type_folder_name = "Labeling"

    output_path = "Root/outputMultiCellulaire/" + study_type_folder_name + "/"

    output_folders_name = [f for f in os.listdir(output_path)]

    bool_diff = ["Yes","No"]
    rn_name = radionuclide_entry.get()
    nb_particles_per_cell = ["1", "2", "3", "4", "5" ,"6", "7", "8", "9" ,"10", "42"]

    type_cell = cell_line_combobox.current()

    available_data_date = []
    available_data_name_file = []
    for i in range(0, len(output_folders_name)):
        if ("_" + cell_compartment + "_" + str(spheroid_compaction) + "CP_" + str(r_sph) + "um_" +
            rn_name + "_diff" + bool_diff[diffusion_combobox.current()] + "_" +
            str(nb_particles_per_cell[number_particles_per_cell_combobox.current()] + "ppc")) in output_folders_name[i]:
            available_data_date.append(output_folders_name[i][0:10])
            available_data_name_file.append(output_folders_name[i])

    available_data_label = tkinter.Label(window, text="Data available :", fg='blue')
    available_data_label.place(x=100, y=510)
    available_data_choice = tkinter.StringVar()
    available_data_combobox = tkinter.ttk.Combobox(window, width=35, textvariable=available_data_choice)
    available_data_combobox['values'] = available_data_date
    available_data_combobox.place(x=400, y=510)

    nom_fichier_root = "output_"
    # Les fichiers root, contenus dans le dossier_root, s'appellent nom_fichier_root{0,...}.root

    validate_button_2 = tkinter.Button(window, text="Validate", command=main)
    validate_button_2.place(x=500, y=560)

    progress_bar = tkinter.ttk.Progressbar(
        window,
        orient='horizontal',
        mode='determinate',
        length=280,
        value=0)
    progress_bar.grid(column=0, row=0, columnspan=2, padx=400, pady=620)

    progress_bar_label = tkinter.ttk.Label(window, text = update_progress_bar_label())
    progress_bar_label.place(x=450, y=650)


if __name__ == '__main__':
    graphic_window()
