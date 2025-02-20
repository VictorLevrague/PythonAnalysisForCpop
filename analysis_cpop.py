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
import nanox_low_energy as nanox
import numpy as np
import os
import pandas as pd
from pathlib import Path
import sys
import uproot
import time
import warnings

warnings.filterwarnings("ignore")
np.set_printoptions(threshold=sys.maxsize)

START_TIME = time.perf_counter()


KEV_IN_J = 1.60218 * 1e-16
WATER_DENSITY = 1e-15  #kg/µm³
UNIT_COEFFICIENT_A = KEV_IN_J / WATER_DENSITY # Gy.µm³.keV-1
ALPHA_PHOTON = [0.313, 0.184, 0.228] #HSG, V79, CHO-K1
BETA_PHOTON = [0.062, 0.020, 0.020]
BETAG = [0.0961, 0.0405, 0.0625]  # constant of Monini et al. 2019

#################### User parameters ####################
CELL_LINE = "V79"
R_SPHEROID = 100 #um
# LIST_PPC = [328, 656, 983, 1311, 1639, 3278, 4917, 6556, 8195]
LIST_PPC = [1639]
#########################################################


def remove_null_values(ei_ef_unique, nb_particles_per_nucleus):
    """
    Removes values = 0 from both arrays
    """
    nb_particles_per_nucleus_np = nb_particles_per_nucleus.to_numpy()
    mask_zero = nb_particles_per_nucleus_np != 0
    return ei_ef_unique.to_numpy()[mask_zero], nb_particles_per_nucleus_np[mask_zero]


def remove_nan_values(df, column):
    """
    Filters out rows where `column` has NaN values.
    """
    return df[pd.notna(df[column])]


def mean_and_std_calculation_dataframe(analysis_dataframe):
    """
    Computes mean and standard deviation for multiple grouped statistics and cleans data.
    """

    columns = [
        'tcp', "tcp_local", "tcp_global", 'cell_survival_local', 'cell_survival_global', 'cell_survival_total',
        'rbe_furusawa', 'rbe_nakano', 'dose_nucleus_(gy)', 'dose_cytoplasm_(gy)', 'dose_cell_(gy)',
        'cross_fire_nucleus', 'nb_particles_per_nucleus', 'edep_spheroid'
    ]

    mean_df = pd.DataFrame()
    mean_df["id_cell"] = analysis_dataframe["id_cell"].unique()

    ei_ef_without_zeros, nb_particles_per_nucleus_without_zeros = remove_null_values(
        analysis_dataframe['ei_ef_sum'], analysis_dataframe['nb_particles_per_nucleus']
    )
    edep_moy_per_nucleus_cross = ei_ef_without_zeros / nb_particles_per_nucleus_without_zeros

    for column in columns:
        mean_df[f"{column}_std"] = analysis_dataframe.groupby(['id_cell'])[column].std()
        mean_df[column] = analysis_dataframe.groupby(['id_cell'])[column].mean()

    mean_df["edep_moy_per_nucleus_cross_(kev)"] = edep_moy_per_nucleus_cross.mean()
    return mean_df


def open_root_file(folder_name, simulation_id):
    """ Opens a ROOT file and extracts relevant simulation data. """

    root_file_name = f"{folder_name}/output_{simulation_id}_t0.root"
    print("root_id =", simulation_id)

    f = uproot.open(root_file_name)

    mandatory_fields = {
        'nameParticle': 'cell/nameParticle',
        'Ei': 'cell/Ei',
        'Ef': 'cell/Ef',
        'ID_Cell': 'cell/ID_Cell',
        'Cellule_D_Emission': 'cell/Cellule_D_Emission',
        'eventID': 'cell/eventID',
        'fEdepn': 'cell/fEdepn',
        'fEdepc': 'cell/fEdepc',
    }
    optional_fields = {
        'indice_if_diffusion': 'cell/indice_if_diffusion',
        'fEdep_sph': 'cell/fEdep_sph',
    }

    # Retrieve required fields
    root_data = {key: f[path].array(library="np") for key, path in mandatory_fields.items()}

    # Retrieve optional fields, handling missing keys
    available_info = {}
    for key, path in optional_fields.items():
        try:
            root_data[key] = f[path].array(library="np")
            available_info[key] = 1
        except uproot.KeyInFileError:
            warnings.warn(f"{key} not available in the dataset.")
            available_info[key] = 0

    root_data_opened = np.core.records.fromarrays(
        list(root_data.values()), names=', '.join(root_data.keys())
    )
    return root_data_opened, available_info.get('indice_if_diffusion', 0), available_info.get('fEdep_sph', 0)


def extract_particle_data(root_data_opened):
    """Extracts event-level and run-level data based on particle types."""
    alpha_particles_mask = np.isin(root_data_opened["nameParticle"], ['alpha', 'alpha+', 'helium'])
    event_data = root_data_opened[alpha_particles_mask]
    run_data_mask = root_data_opened["nameParticle"] == 'EndOfRun'
    run_data = root_data_opened[run_data_mask]
    return event_data, run_data


def modify_ids(event_data, real_id_cells):
    """Maps ID_Cell and Cellule_D_Emission to the real_id_cells."""
    event_data['ID_Cell'] = np.searchsorted(real_id_cells, event_data['ID_Cell'])
    event_data['Cellule_D_Emission'] = np.searchsorted(real_id_cells, event_data['Cellule_D_Emission'])
    return event_data


def compute_energy_loss(event_data, cell_line):
    """Computes energy loss per unique cell."""
    ei, ef = event_data["Ei"], event_data["Ef"]
    id_cells = event_data["ID_Cell"]

    dn1_de_correction = nanox.dn1_de_continuous_mv_tables_global_events_correction(cell_line, "em", "helium",
                                                                                   let="GEANT4",
                                                                                   method_threshold="Interp")

    n1 = nanox.number_of_lethal_events_for_alpha_traversals(dn1_de_correction, np.max(ei))

    energy_loss = n1(ei) - n1(ef)

    # Aggregate per unique cell
    unique_cells = np.unique(id_cells)
    energy_loss_per_cell = np.array([energy_loss[id_cells == cell].sum() for cell in unique_cells])

    return energy_loss_per_cell, nanox.z_tilde_func(ei, ef, cell_line, "helium")


def compute_cross_fire(event_data, available_diffusion):
    ind_non_cross_fire = event_data["ID_Cell"] == event_data["Cellule_D_Emission"]

    if available_diffusion:
        ind_non_cross_fire = ((event_data["ID_Cell"] == event_data["Cellule_D_Emission"]) &
                              (event_data["indice_if_diffusion"] == 0))

    data_noyau_non_cross_fire = event_data[ind_non_cross_fire]
    dose_noyau_non_cross_fire = data_noyau_non_cross_fire["Ei"] - data_noyau_non_cross_fire["Ef"]
    dose_noyau_cross_fire = event_data["Ei"] - event_data["Ef"]
    dose_noyau_cross_fire = np.setdiff1d(dose_noyau_cross_fire, dose_noyau_non_cross_fire)
    sum_dose_noyau_crossfire = np.sum(dose_noyau_cross_fire)
    sum_dose_noyau_non_cross_fire = np.sum(dose_noyau_non_cross_fire)

    sum_dose_noyau_tot = sum_dose_noyau_crossfire + sum_dose_noyau_non_cross_fire
    return sum_dose_noyau_crossfire / sum_dose_noyau_tot


def compute_biological_dose(mean_survival, alpha, beta):
    biological_dose = (np.sqrt(alpha ** 2 - 4 * beta * np.log(np.mean(mean_survival))) - alpha) / (2 * beta)
    return biological_dose


def compute_rbe(mean_survival, spheroid_dose):
    alpha_ref_hsg_aoki_nakano = 0.259
    beta_ref_hsg_aoki_nakano = 0.040
    cell_id = 0
    match CELL_LINE:
        case "HSG":
            cell_id = 0
        case "V79":
            cell_id = 1
        case "CHO-K1":
            cell_id = 2
    dose_bio_furusawa = compute_biological_dose(mean_survival, ALPHA_PHOTON[cell_id], BETA_PHOTON[cell_id])
    print(dose_bio_furusawa)
    dose_bio_nakano = compute_biological_dose(mean_survival, alpha_ref_hsg_aoki_nakano, beta_ref_hsg_aoki_nakano)
    rbe_furusawa = dose_bio_furusawa / spheroid_dose
    print(spheroid_dose)
    print(rbe_furusawa)
    print()
    rbe_nakano = dose_bio_nakano / spheroid_dose
    return rbe_furusawa, rbe_nakano


def compute_doses(run_data, cell_masses):
    """Calculates nuclear, cytoplasmic, and total cell doses."""
    masses_cytoplasms, masses_nuclei, _ = cell_masses
    dosen = run_data["fEdepn"] * KEV_IN_J / masses_nuclei
    dosec = run_data["fEdepc"] * KEV_IN_J / masses_cytoplasms
    return dosen, dosec, dosen + dosec


def compute_survival(n_unique, n_sub_unique):
    """Computes survival probabilities."""
    survival_local = np.exp(-n_unique)
    survival_global = np.exp(-n_sub_unique)
    survival_total = survival_local * survival_global

    # Avoid numerical underflows
    survival_local[survival_local == 0] = 1e-299
    survival_total[survival_total == 0] = 1e-299
    return survival_local, survival_global, survival_total


def compute_tcp(survival_total):
    """Computes Tumor Control Probability (TCP) using Poisson and binomial models."""
    return np.prod(np.exp(-survival_total)), np.prod(1 - survival_total)


def calculations_from_root_file(analysis_dataframe, root_data_opened, real_id_cells, test_file_not_empty,
                                cell_line, masses_file, elements_to_remove, available_diffusion, available_edep_sph):
    """
    Opens root file corresponding to a MC simulation and calculates quantities like cell survival.
    Returns updated analysis dataframe.
    """

    analysis_dataframe_temp = pd.DataFrame()
    event_data, run_data = extract_particle_data(root_data_opened)
    event_data = modify_ids(event_data, real_id_cells)
    analysis_dataframe_temp['id_cell'] = np.arange(len(real_id_cells))

    # Remove redundant cell ids if needed
    if test_file_not_empty:
        run_data = np.delete(run_data, elements_to_remove, 0)

    n_unique, n_sub_tab = compute_energy_loss(event_data, cell_line)
    n_sub_unique = np.bincount(event_data["ID_Cell"].astype(int), weights=n_sub_tab)

    masses_cytoplasms, masses_nuclei, masses_cells = geometry_informations.masses_cells_reading(masses_file)
    dosen, dosec, dose_total = compute_doses(run_data, (masses_cytoplasms, masses_nuclei, masses_cells))
    survival_local, survival_global, survival_total = compute_survival(n_unique, n_sub_unique)
    tcp_poisson, tcp_binomial = compute_tcp(survival_total)
    tcp_poisson_local, tcp_binomial_local = compute_tcp(survival_local)
    tcp_poisson_global, tcp_binomial_global = compute_tcp(survival_global)
    cross_fire_nucleus = compute_cross_fire(event_data, available_diffusion)

    rbe_furusawa, rbe_nakano = 0, 0
    r_tum = float(R_SPHEROID) * 10**(-6) #in meters
    masse_tum = ((4/3)*np.pi*r_tum**3)*1000 #in kg
    spheroid_dose = run_data[0]["fEdep_sph"] * KEV_IN_J / masse_tum
    if available_edep_sph:
        rbe_furusawa, rbe_nakano = compute_rbe(survival_total, spheroid_dose)

    nb_particles_per_nucleus = np.bincount((event_data["ID_Cell"]).astype(int))
    ei_ef_unique_sur_une_simu = np.bincount(event_data["ID_Cell"].astype(int),
                                            weights=event_data["Ei"] - event_data["Ef"])

    # Store results
    analysis_dataframe_temp['dose_nucleus_(gy)'] = dosen
    analysis_dataframe_temp['dose_cytoplasm_(gy)'] = dosec
    analysis_dataframe_temp['dose_cell_(gy)'] = dose_total
    analysis_dataframe_temp['cell_survival_local'] = survival_local
    analysis_dataframe_temp['cell_survival_global'] = survival_global
    analysis_dataframe_temp['cell_survival_total'] = survival_total
    analysis_dataframe_temp['tcp'] = tcp_binomial
    analysis_dataframe_temp['tcp_local'] = tcp_binomial_local
    analysis_dataframe_temp['tcp_global'] = tcp_binomial_global
    analysis_dataframe_temp['rbe_furusawa'] = rbe_furusawa
    analysis_dataframe_temp['rbe_nakano'] = rbe_nakano
    analysis_dataframe_temp['cross_fire_nucleus'] = cross_fire_nucleus
    analysis_dataframe_temp['ei_ef_sum'] = ei_ef_unique_sur_une_simu
    analysis_dataframe_temp['nb_particles_per_nucleus'] = nb_particles_per_nucleus
    analysis_dataframe_temp['edep_spheroid'] = spheroid_dose

    return pd.concat([analysis_dataframe, analysis_dataframe_temp], ignore_index=True)


def create_folder_for_output_analysis_files(folder_name):
    base_path = os.path.join(folder_name)
    try:
        os.makedirs(base_path, exist_ok=True)
        print(f"Successfully created or verified: {base_path}")
    except Exception as error:
        print(f"Error creating folder {base_path}: {error}")


def analysis(folder_root, folder_analysis, xml_geom, name_geom, nb_simulations, masses_file):
    """Main function to analyze ROOT files and generate CSV reports."""

    start_time = time.perf_counter()
    print("folder_root:", folder_root)
    print("folder_analysis:", folder_analysis)
    create_folder_for_output_analysis_files(folder_analysis)

    nb_cells_xml = geometry_informations.count_number_of_cells_in_xml_file(xml_geom)

    txt_id_deleted_cells = Path(f"Cpop_Deleted_Cells_ID_Txt/Previous_Data/IDCell_{name_geom}.txt")
    real_id_cells, test_file_not_empty, deleted_id_txt = geometry_informations.cpop_real_cell_id_determination(
        txt_id_deleted_cells, nb_cells_xml)

    root_data, available_diffusion, available_edep_sph = open_root_file(folder_root, 0)
    elements_to_remove = remove_cell_ids(root_data, test_file_not_empty, deleted_id_txt, real_id_cells)

    analysis_dataframe = pd.DataFrame()
    for sim_id in range(nb_simulations):
        analysis_dataframe = process_root_file(folder_root, sim_id, analysis_dataframe, real_id_cells,
                                               test_file_not_empty, elements_to_remove, masses_file,
                                               available_diffusion, available_edep_sph)
    save_analysis_results(analysis_dataframe, folder_analysis)
    log_execution_time(start_time)


def process_root_file(folder_name, simulation_id, analysis_dataframe, real_id_cells, test_file_not_empty,
                      elements_to_remove, masses_file, available_diffusion, available_edep_sph):
    """Processes a single ROOT file and updates the analysis DataFrame."""

    root_data, _, _ = open_root_file(folder_name, simulation_id)
    analysis_dataframe = calculations_from_root_file(analysis_dataframe, root_data, real_id_cells,
                                                     test_file_not_empty, CELL_LINE, masses_file, elements_to_remove,
                                                     available_diffusion, available_edep_sph)
    return analysis_dataframe


def remove_cell_ids(root_data_opened, test_file_not_empty, deleted_id_txt, real_id_cells):
    mapped_cell_ids = np.searchsorted(real_id_cells, np.unique(real_id_cells))

    event_particles = np.isin(root_data_opened["nameParticle"], ['alpha', 'alpha+', 'helium'])
    data_event_level = root_data_opened[event_particles]

    end_of_run_mask = root_data_opened["nameParticle"] == 'EndOfRun'
    data_run_level = root_data_opened[end_of_run_mask]

    # Map ID_Cell and Cellule_D_Emission to corrected IDs
    data_event_level = data_event_level.copy()
    data_event_level["ID_Cell"] = mapped_cell_ids[np.searchsorted(real_id_cells, data_event_level["ID_Cell"])]
    data_event_level["Cellule_D_Emission"] = mapped_cell_ids[
        np.searchsorted(real_id_cells, data_event_level["Cellule_D_Emission"])]

    elements_to_remove = []
    if test_file_not_empty:
        elements_to_remove = np.where(np.isin(data_run_level["ID_Cell"], deleted_id_txt))[0].tolist()

    return elements_to_remove


def save_analysis_results(analysis_dataframe, folder_analysis):
    """Saves the processed analysis data to CSV files."""
    results_folder = Path(f"{folder_analysis}")
    results_folder.mkdir(parents=True, exist_ok=True)

    analysis_dataframe.to_csv(results_folder / "AllData.csv", index=False)
    mean_and_std_calculation_dataframe(analysis_dataframe).to_csv(results_folder / "Results.csv", index=False)


def log_execution_time(start_time):
    """Logs the total execution time of the script."""
    end_time = time.perf_counter()
    total_seconds = int(end_time - start_time)
    print(f"\nTotal time: {total_seconds // 60} minutes {total_seconds % 60} seconds")


def main():
    study_name = "NetiComparison"
    geometry_name = "Net100um47CP"
    distribution_type = "Uniform"
    xml_geometry_file = f"Cpop_Geom_XML/{geometry_name}.cfg.xml"
    radionuclide = "Po210"
    labeling = 1
    date = "2025_02_19"
    nb_simulations = 20
    ####################
    name_config = f"{study_name}/{geometry_name}/{distribution_type}/Labeling{labeling}%/{radionuclide}/{date}"
    analysis_config = f"AnalysisResults/{name_config}"
    masses_file = f"Cpop_Masse_Txt/MassesCell_{geometry_name}.txt"
    ####################
    for ppc in LIST_PPC:
        print(ppc, " ppc", " + ", distribution_type, " distribution")
        folder_root = f"Root/output/{name_config}/output_{ppc}ppc"
        folder_analysis = f"{analysis_config}/output_{ppc}ppc"
        analysis(folder_root, folder_analysis, xml_geometry_file, geometry_name, nb_simulations, masses_file)


if __name__ == "__main__":
    main()
