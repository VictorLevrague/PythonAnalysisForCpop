import numpy as np
import os
from xml.dom import minidom

def column(matrix, i):
    return [row[i] for row in matrix]

def subset_sorted_array(A,B):
    """
    Removes elements of A sorted array that are contained in B
    """
    Aa = A[np.where(A <= np.max(B))]
    Bb = (B[np.searchsorted(B,Aa)] !=  Aa)
    Bb = np.pad(Bb,(0,A.shape[0]-Aa.shape[0]), mode='constant', constant_values=True)
    return A[Bb]

def count_number_of_cells_in_xml_file(xml_filename):
    xml_file_opened = minidom.parse(xml_filename)
    xml_file_opened_with_cell_tag = xml_file_opened.getElementsByTagName('CELL')
    nb_cells = 0
    for node in xml_file_opened_with_cell_tag:
        xml_file_opened_with_cell_and_x_tag = node.getElementsByTagName('x')
        for xml_parser in xml_file_opened_with_cell_and_x_tag:
            # x = xml_parser.firstChild.data
            nb_cells += 1
    return(nb_cells)

def cpop_real_cell_id_determination(file_with_deleted_cell_id_in_cpop, nb_cellules_xml):
    """
    When geometry is generated in CPOP, cell ids are between 3 and nb_cells+3.
    But, cells are removed during the cell overlap mangement,
    hence, the associated cell ids don't exist anymore.

    Returns
    ------
    a sorted array with the cell ids that exist in the geometry
    a test for the validity of the txt file for deleted ids
    a sorted array with the deleted cells ids
    """
    test_file_not_empty = os.stat(file_with_deleted_cell_id_in_cpop).st_size
    if test_file_not_empty != 0:
        deleted_id_txt = np.loadtxt(file_with_deleted_cell_id_in_cpop)
        deleted_id_txt = np.unique(deleted_id_txt)
        deleted_id_txt = np.sort(deleted_id_txt)
    real_id_cells = np.arange(3,nb_cellules_xml+3)
    if test_file_not_empty != 0:
        real_id_cells = subset_sorted_array(real_id_cells,deleted_id_txt)
    return(real_id_cells, test_file_not_empty, deleted_id_txt)

def masses_cells_reading(txt_file_with_masses_cells):
    """
    Returns 3 numpy arrays with masses of nuclei, cytoplams and cells, sorted by cell ids
    """
    masses_cells_txt_numpy = np.loadtxt(txt_file_with_masses_cells, dtype={'names': ('masse_noyau', 'unit1', 'masse_cell', 'unit2'),
                                                    'formats': (float, '|S15', float, '|S15')})
    masses_nuclei = (column(masses_cells_txt_numpy,0))
    masses_nuclei = np.array(masses_nuclei) * 10**(-6) #conversion in kg
    masses_cells = (column(masses_cells_txt_numpy,2))
    masses_cells = np.array(masses_cells) * 10 ** (-6)  #conversion in kg
    masses_cytoplasms = masses_cells - masses_nuclei  #kg
    return(masses_cytoplasms, masses_nuclei, masses_cells)

def positions_cells_reading(xml_file_with_cells_positions, real_id_cells):
    """
    Returns 3 numpy arrays, with the x, y and z positions of the cells, sorted by cell ids
    """
    cells_positions_xml_opened = minidom.parse(xml_file_with_cells_positions)
    cells_positions_xml_opened_with_cell_tag = cells_positions_xml_opened.getElementsByTagName('CELL')
    nb_cellules_xml = count_number_of_cells_in_xml_file(xml_file_with_cells_positions)
    ###
    positions_x=np.zeros(nb_cellules_xml)
    positions_y=np.zeros(nb_cellules_xml)
    positions_z=np.zeros(nb_cellules_xml)
    positions_and_id = np.zeros((nb_cellules_xml,4))
    row_nb = 0
    ###
    for parser_xml in cells_positions_xml_opened_with_cell_tag:
        positions_and_id[row_nb][3] = parser_xml.attributes['ID'].value
        row_nb += 1
    row_nb = 0
    for node in cells_positions_xml_opened_with_cell_tag:
        positions_x_xml = node.getElementsByTagName('x')
        positions_y_xml = node.getElementsByTagName('y')
        positions_z_xml = node.getElementsByTagName('z')
        for parser_xml in positions_x_xml:
            positions_and_id[row_nb][0] = parser_xml.firstChild.data
        for parser_xml in positions_y_xml:
            positions_and_id[row_nb][1] = parser_xml.firstChild.data
        for parser_xml in positions_z_xml:
            positions_and_id[row_nb][2] = parser_xml.firstChild.data
        row_nb += 1
    positions_and_id = positions_and_id[positions_and_id[:, 3].argsort()] #sorts the array by cells id
    index_cell_in_positions_and_id = 0
    indexes_to_delete = []
    for cell_id in column(positions_and_id,3):
        if not ((cell_id in real_id_cells)):
            indexes_to_delete.append(index_cell_in_positions_and_id)
        index_cell_in_positions_and_id += 1
    positions_and_id = np.delete(positions_and_id, indexes_to_delete, 0)
    positions_x = positions_and_id[:,0]
    positions_y = positions_and_id[:,1]
    positions_z = positions_and_id[:,2]
    return(positions_x, positions_y, positions_z)

geometry_name = "Elg095um75CP"
###
txt_cells_masses = "Cpop_Masse_Txt/" + "MassesCell_" + geometry_name + ".txt"
try :
    os.makedirs(os.path.join("./GeometryInformations/" + geometry_name))
except:
    pass
masses_cytoplasms, masses_nuclei, masses_cells = masses_cells_reading(txt_cells_masses)
np.savetxt("GeometryInformations/" + geometry_name + "/MassesCells.txt", (masses_cytoplasms,
                                                                 masses_nuclei, masses_cells))
###
xml_geometry_file = "Cpop_Geom_XML/" + geometry_name + ".cfg" + ".xml"
nb_cellules_xml = count_number_of_cells_in_xml_file(xml_geometry_file)
txt_id_deleted_cells = "Cpop_Deleted_Cells_ID_Txt/" + "IDCell_" + geometry_name + ".txt"
real_id_cells, test_file_not_empty, deleted_id_txt = cpop_real_cell_id_determination(txt_id_deleted_cells,
                                                                                     nb_cellules_xml)
positions_x, positions_y, positions_z = positions_cells_reading(xml_geometry_file, real_id_cells)
np.savetxt("GeometryInformations/" + geometry_name + "/PositionsCells.txt", (positions_x,
                                                                 positions_y, positions_z))