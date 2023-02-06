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

############### Fit parameters used in article TCP RIV-alpha ###################
Y0 = [0.06072969, 0.02562553, 0.03934994] #HSG, V79, CHO-K1
A = [-0.18385472, -0.10426184, -0.11163773]
W = [3.05093045, 2.87758559, 3.20398251]
XC = [0.46545609, 0.38084839, 0.48452192]

# ################ Fit parameters from 2022/12/16 ################# TO DO : evaluate differences with previous parameters
# Y0 = [0.06486164, 0.02722410, 0.04221387] #HSG, V79, CHO-K1
# A = [-0.26336407, -0.11801719, -0.19357751]
# W = [3.39940424, 2.97713123, 3.90866411]
# XC = [-0.00863166, 0.23348883, -0.25238105]

## TO DO : generalize alpha & beta photon for all available cell lines (HSG, V79 & CHO-K1)
ALPHA_PHOTON = 0.313  # HSG
BETA_PHOTON = 0.0615  # HSG

BETAG = [0.0961, 0.0405, 0.0625]  # constant of Monini et al. 2019


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

    tables_conversion_energy_in_let = pd.read_excel(f"E_TEL/conversion_tables_{data_base}.xlsx").to_records()
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

def deletion_of_cells_with_no_alpha_traversals(ei_ef_unique, nb_particles_per_nucleus):
    """
    Deletion of cells that have not been crossed by alpha particles,
    for the calculation of the mean deposited energy per nucleus cross
    TO DO : complete docstring and rewrite this function
    """

    ei_ef_unique_np = ei_ef_unique.to_numpy()
    nb_particles_per_nucleus_np = nb_particles_per_nucleus.to_numpy()

    ind_nucleus_with_no_particles = [index for index, value in enumerate(
        nb_particles_per_nucleus_np) if value == 0]

    nb_particles_per_nucleus_without_zeros = np.array([])
    ei_ef_without_zeros = np.array([])

    nb_particles_per_nucleus_without_zeros =\
        np.append(nb_particles_per_nucleus_without_zeros,
                  np.delete( nb_particles_per_nucleus_np, ind_nucleus_with_no_particles))
    ei_ef_without_zeros = np.append(ei_ef_without_zeros,
                  np.delete(ei_ef_unique_np, ind_nucleus_with_no_particles))
    return ei_ef_without_zeros, nb_particles_per_nucleus_without_zeros

def mean_and_std_calculation_dataframe(analysis_dataframe):
    """
    TO DO : docstring
    """

    ei_ef_without_zeros, nb_particles_per_nucleus_without_zeros = deletion_of_cells_with_no_alpha_traversals(
                analysis_dataframe['ei_ef_sum'], analysis_dataframe['nb_particles_per_nucleus'])


    edep_moy_per_nucleus_cross = ei_ef_without_zeros / nb_particles_per_nucleus_without_zeros
    analysis_dataframe['tcp'] = analysis_dataframe.groupby(['id_cell'])['tcp'].mean()
    analysis_dataframe['cell_survival_local_std'] = analysis_dataframe.groupby(['id_cell'])\
                                                                        ['cell_survival_local'].std()
    analysis_dataframe['cell_survival_local'] = analysis_dataframe.groupby(['id_cell'])['cell_survival_local'].mean()
    analysis_dataframe['cell_survival_global_std'] = analysis_dataframe.groupby(['id_cell'])\
                                                                        ['cell_survival_global'].std()
    analysis_dataframe['cell_survival_global'] = analysis_dataframe.groupby(['id_cell'])['cell_survival_global'].mean()
    analysis_dataframe['dose_nucleus_std'] = analysis_dataframe.groupby(['id_cell']) \
                                                                        ['dose_nucleus_(gy)'].std()
    analysis_dataframe['dose_nucleus_(gy)'] = analysis_dataframe.groupby(['id_cell'])['dose_nucleus_(gy)'].mean()
    analysis_dataframe['dose_cytoplasm_std'] = analysis_dataframe.groupby(['id_cell']) \
                                                                        ['dose_cytoplasm_(gy)'].std()
    analysis_dataframe['dose_cytoplasm_(gy)'] = analysis_dataframe.groupby(['id_cell'])['dose_cytoplasm_(gy)'].mean()
    analysis_dataframe['dose_cell_(gy)'] = analysis_dataframe.groupby(['id_cell'])['dose_cell_(gy)'].mean()
    analysis_dataframe['cross_fire_nucleus'] = analysis_dataframe.groupby(['id_cell'])['cross_fire_nucleus'].mean()
    analysis_dataframe['cross_fire_nucleus_zone1'] = analysis_dataframe.groupby(['id_cell'])\
                                                                        ['cross_fire_nucleus_zone1'].mean()
    analysis_dataframe['cross_fire_nucleus_zone2'] = analysis_dataframe.groupby(['id_cell'])\
                                                                        ['cross_fire_nucleus_zone2'].mean()
    analysis_dataframe['tcp_formula_poisson'] = analysis_dataframe.groupby(['id_cell'])['tcp_formula_poisson'].mean()
    analysis_dataframe['biological_dose_(gy)'] = analysis_dataframe.groupby(['id_cell'])['biological_dose_(gy)'].mean()
    analysis_dataframe['edep_moy_per_nucleus_cross_(kev)'] = edep_moy_per_nucleus_cross.mean()

    analysis_dataframe.drop('ei_ef_sum', inplace=True, axis=1)
    analysis_dataframe.drop('nb_particles_per_nucleus', inplace=True, axis=1)

    return analysis_dataframe

def open_root_file(simulation_id):
    root_file_name = f"Root/outputMultiCellulaire/{dossier_root}{nom_fichier_root}{simulation_id}_t0.root"

    if verbose == 1:
        print("Root file name : ", root_file_name)

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

    return root_data_opened, indice_available_diffusion_info, indice_available_edep_sph_info

def calculations_from_root_file(analysis_dataframe, root_data_opened, simulation_id, indice_available_diffusion_info,
                                indice_available_edep_sph_info, real_id_cells, test_file_not_empty, deleted_id_txt):
    """
    Opens root file corresponding to a MC simulation and calculates quantities like cell survivals
    Returns Check
    """

    analysis_dataframe_temp = pd.DataFrame()

    nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus = []
    nb_cellules_reel = len(real_id_cells)
    perfect_id_cells = np.arange(0,nb_cellules_reel)

    # analysis_dataframe_temp['simulation_id'] = np.full((1, nb_cellules_reel), simulation_id)[0]

    analysis_dataframe_temp['id_cell'] = np.arange(nb_cellules_reel)
    analysis_dataframe_temp['zone_cell'] = zone_cell

    ind_alphaplusplus = root_data_opened["nameParticle"] == 'alpha'
    ind_alphaplus = root_data_opened["nameParticle"] == 'alpha+'
    ind_helium = root_data_opened["nameParticle"] == 'helium'

    data_event_level = (np.concatenate((root_data_opened[ind_alphaplusplus],
                                       root_data_opened[ind_alphaplus],
                                       root_data_opened[ind_helium])))

    ind_end_of_run = root_data_opened["nameParticle"] == 'EndOfRun'

    data_edep_cell = root_data_opened[ind_end_of_run]

    ########################## Vérification diffusion aux bonnes énergies ###############################

    if indice_available_diffusion_info == 1:

        unique_data_event_level_event_id = np.unique(data_event_level['eventID'], return_index=True)

        ind_diff_0 = ind_diff_1 = len_unique = 0

        indices_ab = unique_data_event_level_event_id[1]

        unique_data_event_level_ind_diff_corresponding_to_unique_event_id = \
            np.take(data_event_level['indice_if_diffusion'], indices_ab)

        for i in range(0, len(unique_data_event_level_ind_diff_corresponding_to_unique_event_id)):
            if unique_data_event_level_ind_diff_corresponding_to_unique_event_id[i] == 0:
                ind_diff_0 += 1
                len_unique += 1
            elif unique_data_event_level_ind_diff_corresponding_to_unique_event_id[i] == 1:
                ind_diff_1 += 1
                len_unique += 1

        if verbose == 1:
            print("% d'event sans diffusion = ", ind_diff_0 / len_unique)
            print("% d'event avec diffusion = ", ind_diff_1 / len_unique)

    ####################### Modification des ID de CPOP ###################################

    ################ data_event_level #########################################
    for ind_modif_id in range(0, len(data_event_level)):
        index_id_cell = np.where(real_id_cells == data_event_level[ind_modif_id]["ID_Cell"])
        data_event_level[ind_modif_id]["ID_Cell"] = perfect_id_cells[index_id_cell]

        index_cellule_emission = np.where(real_id_cells == data_event_level[ind_modif_id]["Cellule_D_Emission"])
        data_event_level[ind_modif_id]["Cellule_D_Emission"] = perfect_id_cells[index_cellule_emission]

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

    ei = data_event_level["Ei"]  # Energy in keV
    ef = data_event_level["Ef"]
    dn1_de_continous_pre_calculated = dn1_de_continous()
    n1 = number_of_lethal_events_for_alpha_traversals(dn1_de_continous_pre_calculated)

    n_tab = (n1(ei) - n1(ef))

    txt_cells_masses=f"Cpop_Masse_Txt/MassesCell_{nom_config}.txt"
    masses_cytoplasms, masses_nuclei, masses_cells = geometry_informations.masses_cells_reading(txt_cells_masses)
    dosen_append_sur_une_simu_np = (((data_edep_cell)["fEdepn"]) * KEV_IN_J / masses_nuclei)
    dosec_append_sur_une_simu_np = (((data_edep_cell)["fEdepc"]) * KEV_IN_J / masses_cytoplasms)

    analysis_dataframe_temp['dose_nucleus_(gy)'] = dosen_append_sur_une_simu_np
    analysis_dataframe_temp['dose_cytoplasm_(gy)'] = dosec_append_sur_une_simu_np
    analysis_dataframe_temp['dose_cell_(gy)'] = dosen_append_sur_une_simu_np + dosec_append_sur_une_simu_np

    nb_particles_per_nucleus = np.bincount((data_event_level["ID_Cell"]).astype(int))
    ei_ef_unique_sur_une_simu  = np.bincount(data_event_level["ID_Cell"].astype(int), weights=ei - ef)

    while len(ei_ef_unique_sur_une_simu) < nb_cellules_reel:
        ei_ef_unique_sur_une_simu = np.append(ei_ef_unique_sur_une_simu, 0)

    while len(nb_particles_per_nucleus) < nb_cellules_reel:
        nb_particles_per_nucleus = np.append(nb_particles_per_nucleus, 0)

    #################################### Cross-fire dose au noyau ##########################################

    ############################################

    ind_non_cross_fire = data_event_level["ID_Cell"] == data_event_level["Cellule_D_Emission"]
    ind_cross_fire = data_event_level["ID_Cell"] != data_event_level["Cellule_D_Emission"]

    if indice_available_diffusion_info == 1:
        ind_non_cross_fire = ((data_event_level["ID_Cell"] == data_event_level["Cellule_D_Emission"]) &
                              (data_event_level["indice_if_diffusion"] == 0))
        ind_cross_fire = ((data_event_level["ID_Cell"] != data_event_level["Cellule_D_Emission"]) &
                          (data_event_level["indice_if_diffusion"] == 1))

    data_noyau_non_cross_fire = data_event_level[ind_non_cross_fire]

    dose_noyau_non_cross_fire = data_noyau_non_cross_fire["Ei"] - data_noyau_non_cross_fire["Ef"]

    dose_noyau_cross_fire = data_event_level["Ei"] - data_event_level["Ef"]

    dose_noyau_cross_fire = np.setdiff1d(dose_noyau_cross_fire, dose_noyau_non_cross_fire)

    sum_dose_noyau_crossfire = np.sum(dose_noyau_cross_fire)
    sum_dose_noyau_non_cross_fire = np.sum(dose_noyau_non_cross_fire)

    ################################################
    indice_zone1 = zone_cell == 1
    indice_zone2 = zone_cell == 2

    ei_ef_unique_non_cross_fire = np.bincount(((data_event_level["ID_Cell"])[ind_non_cross_fire]).astype(int),
                                              weights=ei[ind_non_cross_fire] - ef[ind_non_cross_fire])
    ei_ef_unique_cross_fire = np.bincount(((data_event_level["ID_Cell"])[ind_cross_fire]).astype(int),
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

    count_event_id = np.bincount((data_event_level["eventID"]).astype(int))

    for id_part in range(0, len(np.unique(data_event_level["eventID"]))):
        nb_nucl_traversees_par_la_particule = count_event_id[id_part]
        nb_nucl_traversees_par_la_particule_tab_sur_toutes_simus.append(nb_nucl_traversees_par_la_particule)

    # print("############################# Calcul de survie #######################################")

    n_unique = np.bincount(data_event_level["ID_Cell"].astype(int), weights=n_tab)
    while len(n_unique) < nb_cellules_reel:
        n_unique = np.append(n_unique, 0)

    n_unique_tot_sur_une_simu = n_unique

    progress_bar['value'] += round(100 / (nb_complete_simulations - nb_files_with_errors), 2)
    progress_bar_label['text'] = update_progress_bar_label()
    window.update_idletasks()

    sum_dose_noyau_tot = sum_dose_noyau_crossfire + sum_dose_noyau_non_cross_fire
    analysis_dataframe_temp['cross_fire_nucleus'] = sum_dose_noyau_crossfire / sum_dose_noyau_tot
    sum_dose_noyau_tot_zone1 = sum_dose_noyau_crossfire_zone1 + sum_dose_noyau_non_cross_fire_zone1
    analysis_dataframe_temp['cross_fire_nucleus_zone1'] = sum_dose_noyau_crossfire_zone1 / sum_dose_noyau_tot_zone1
    sum_dose_noyau_tot_zone2 = sum_dose_noyau_crossfire_zone2 + sum_dose_noyau_non_cross_fire_zone2
    analysis_dataframe_temp['cross_fire_nucleus_zone2'] = sum_dose_noyau_crossfire_zone2 / sum_dose_noyau_tot_zone2

    analysis_dataframe_temp['ei_ef_sum'] = ei_ef_unique_sur_une_simu
    analysis_dataframe_temp['nb_particles_per_nucleus'] = nb_particles_per_nucleus

    surviel_append_sur_une_simu = np.exp(-n_unique_tot_sur_une_simu)

    surviel_append_sur_une_simu[np.where(surviel_append_sur_une_simu == 0)] = 10 ** (-299)

    analysis_dataframe_temp['cell_survival_local'] = surviel_append_sur_une_simu

    dose_bio_append_sur_une_simu = \
       (np.sqrt(ALPHA_PHOTON ** 2 - 4 * BETA_PHOTON * np.log(surviel_append_sur_une_simu)) - ALPHA_PHOTON) \
            / (2 * BETA_PHOTON)

    analysis_dataframe_temp['biological_dose_(gy)'] = dose_bio_append_sur_une_simu

    exp_surviel = np.exp(-np.asarray(surviel_append_sur_une_simu))
    tcp_une_simu = np.prod(exp_surviel)
    tcp_test_formula = np.prod(1 - surviel_append_sur_une_simu)
    analysis_dataframe_temp['tcp'] = tcp_une_simu
    analysis_dataframe_temp['tcp_formula_poisson'] = tcp_test_formula

    survieg_append_sur_une_simu = \
        np.exp(-n_unique_tot_sur_une_simu - BETAG[type_cell] * (dosen_append_sur_une_simu_np ** 2))

    analysis_dataframe_temp['cell_survival_global'] = survieg_append_sur_une_simu

    spheroid_dose = data_edep_cell[0]["fEdep_sph"] * KEV_IN_J / masse_tum
    analysis_dataframe_temp['spheroid_dose'] = spheroid_dose


    return pd.concat([analysis_dataframe, analysis_dataframe_temp], ignore_index=True)

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

def print_mean_results(analysis_dataframe):
    print("Mean values with standard deviations : ")
    for name, values in analysis_dataframe.iteritems():
        if name not in ["simulation_id", "id_cell", "zone_cell"]:
            name = name.replace("_", " ").capitalize()
            mean_value = values.mean()
            print(f"{name} : {mean_value}")
    return None

def main():
    global nb_cellules_reel, masses_cells, masse_tum, nb_cell_zone_1, nb_cell_zone_2, nb_files_with_errors, zone_cell,\
    verbose
    verbose = verbose_radiovalue.get()
    create_folder_for_output_analysis_files()
    ##################### Gestion des ID de CPOP ##################################################
    txt_id_deleted_cells = f"Cpop_Deleted_Cells_ID_Txt/IDCell_{nom_config}.txt"
    real_id_cells, test_file_not_empty, deleted_id_txt =\
        geometry_informations.cpop_real_cell_id_determination(txt_id_deleted_cells, nb_cellules_xml)
    nb_cellules_reel = len(real_id_cells)
    ###################### Lecture Geométrie #####################################
    ###### Masses #######
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

    indexes_root_files_without_errors_np, nb_files_with_errors = id_deletion_of_root_outputs_with_errors()

    ##################################################

    analysis_dataframe = pd.DataFrame()

    for simulation_id in indexes_root_files_without_errors_np:
        root_data_np, indice_available_diffusion_info, indice_available_edep_sph_info = open_root_file(simulation_id)
        analysis_dataframe = calculations_from_root_file(analysis_dataframe, root_data_np, simulation_id,
                                                       indice_available_diffusion_info, indice_available_edep_sph_info,
                                                       real_id_cells, test_file_not_empty, deleted_id_txt)

    progress_bar['value'] = math.floor(progress_bar['value'])

    mean_and_std_calculation_dataframe(analysis_dataframe).to_csv('Test_df_to_pandas.csv')

    print()
    print_mean_results(analysis_dataframe)

    # if indice_available_edep_sph_info:
    #     print("Dose sphéroïde")
        # print(spheroid_dose[0], "Gy") #Check : add spheroid dose

    print("Nombre de simus fonctionnelles = ")
    print(nb_complete_simulations - nb_files_with_errors)

    print()

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