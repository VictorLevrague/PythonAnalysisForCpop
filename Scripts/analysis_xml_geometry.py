"""
Input :
- MassesCells.Txt that contains masse of all cells in a CPOP simulation
- Radius of the spheroid

Returns :
The % of cell packing
"""

import geometry_informations
import numpy as np

nom_config = "Grid50CP"
r_sph = 100 #um
folder_name = "Previous_Data"
PATH_TO_ROOT_CPOP_ANALYSIS_FOLDER = "/home/levrague/Documents/Python/Python_these/CPOP/Analyse_CPOP"

def compute_radius_from_mass(masses):
    rho = 10**(-15) #kg/µm3
    radii = ((3*masses)/(rho*4*np.pi))**(1/3)
    return radii

def compute_cell_packing(masses_cells):
    r_tum = float(r_sph) * 10 ** (-6)  # in meters
    masse_tum = ((4 / 3) * np.pi * r_tum ** 3) * 1000  # in kg
    return sum(masses_cells)/masse_tum

def compute_cell_deformation_factor(masses_cells, maximum_radius):
    maximum_masse = ((4 / 3) * np.pi * maximum_radius ** 3) * 1000 #kg
    mean_deformed_masse = np.mean(masses_cells) #kg
    return 1 - mean_deformed_masse / maximum_masse

def main():
    txt_cells_masses = f"{PATH_TO_ROOT_CPOP_ANALYSIS_FOLDER}/Cpop_Masse_Txt/" + folder_name + f"/MassesCell_{nom_config}.txt"
    masses_cytoplasms, masses_nuclei, masses_cells = geometry_informations.masses_cells_reading(txt_cells_masses) #kg

    radii = compute_radius_from_mass(masses_nuclei) #µm
    print("mean nucleus radius: ", np.mean(radii), " += ", np.std(radii) ," µm")

    cell_packing = compute_cell_packing(masses_cells)

    print("Cell packing = ", cell_packing *100, "%")

    maximum_cell_radius = 6 * 10**(-6) #m
    cell_deformation_factor = compute_cell_deformation_factor(masses_cells, maximum_cell_radius)
    print("Cell deformation factor = ", cell_deformation_factor*100, "%")

    maximum_nucleus_radius = 4 * 10**(-6) #m
    nucleus_deformation_factor = compute_cell_deformation_factor(masses_nuclei, maximum_nucleus_radius)
    print("Nucleus deformation factor = ", nucleus_deformation_factor*100, "%")


if __name__ == "__main__":
    main()