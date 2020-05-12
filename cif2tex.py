from pymatgen.io.cif import CifParser
from pymatgen import Structur
import sys, getopt
import pandas as pd

def cif_import(filename):
  parser = CifParser(filename + ".cif")
  structure = parser.as_dict()
  structure = structure[list(structure.keys())[0]]
  return structure

def ciftable_tex(filename, structure):
  # Open the file with writing permission
  with open(filename +".tex", "w") as myfile:

    # Write a line to the file
    myfile.write("\\begin{longtable}\n \\caption{}\n \label{tab:cryst_data} \n")
    myfile.write("\\toprule\n & \\ce{}                               \\* \\midrule \n \\endfirstheadl \n \\endhead \n % \n \\bottomrule \n \\endfoot \n % \n\\endlastfoot \n% \n")
    myfile.write("Formula weight & "+ structure["_chemical_formula_weight"]                               +"    \\ \n")
    myfile.write("Temperature/K & "+ structure["_cell_measurement_temperature"] +"                                    \\ \n")
    myfile.write("Crystal system  & "+ structure["_space_group_crystal_system"]     +"                           \\ \n")
    myfile.write("Space group & "+ structure["_space_group_name_H-M_alt"].replace(" ", "")+ "                                   \\ \n")
    myfile.write("a/\\AA\\ & "+ structure["_cell_length_a"] +"                                \\ \n")
    myfile.write("b/\\AA\\ & "+ structure["_cell_length_b"] +"                              \\ \n")
    myfile.write("c/\\AA\\ & "+ structure["_cell_length_c"] +"                                \\ \n")
    myfile.write("$\\alpha$/\\textdegree  & "+ structure["_cell_angle_alpha"] +"                                        \\ \n")
    myfile.write("$\\beta$/\\textdegree  & "+ structure["_cell_angle_beta"] +"                               \\ \n")
    myfile.write("$\\gamma$/\\textdegree  & "+ structure["_cell_angle_gamma"] +"                                        \\ \n")
    myfile.write("Volume/\\AA$^3$   & "+ structure["_cell_volume"] +"                               \\ \n")
    myfile.write("Z & "+ structure["_cell_formula_units_Z"] +"                                        \\ \n")
    myfile.write("$\\rho$calcg/cm$^3$  & "+ structure["_exptl_crystal_density_meas"]+ "       \\ \n")
    myfile.write("F(000)  & "+ structure["_exptl_crystal_F_000"] +" \\ \n")
    myfile.write("Crystal size/mm$^3$ & "+ structure["_exptl_crystal_size_max"] +" x" + structure["_exptl_crystal_size_mid"] +" x" + structure["_exptl_crystal_size_min"] +"\\\ \n")
    myfile.write("Radiation  & "+ structure["_diffrn_radiation_type"]+ "$\alpha$ ($\\lambda$ = " + structure["_diffrn_radiation_wavelength"] +")                        \\\ \n")
    myfile.write("2$\\theta$ range for data collection/\textdegree                 & "+ structure["_diffrn_reflns_theta_min"]+ " to " + structure["_diffrn_reflns_theta_max"] +" \\\ \n")
    myfile.write("Index ranges & " +structure["_diffrn_reflns_limit_h_min"]+" $\leq$ h $\leq$ "+structure["_diffrn_reflns_limit_h_max"]+","+structure["_diffrn_reflns_limit_k_min"]+"  $\leq$ k $\leq$"+ structure["_diffrn_reflns_limit_k_max"]+"," + structure["_diffrn_reflns_limit_l_min"]+" $\leq$ l $\leq$"+ structure["_diffrn_reflns_limit_h_min"]+"  \\ \n")
    myfile.write("Reflections collected & "+ structure["_diffrn_reflns_number"]+" \\ \n")
    myfile.write("Independent reflections & " + structure["_reflns_number_total"] + "{[}R$_{int}$ =" + structure["_diffrn_reflns_av_R_equivalents"]+ ", R$_\sigma$ =" + structure["_diffrn_reflns_av_unetI/netI"]+"{]} \\ \n")
    myfile.write("Data/restraints/parameters &"+ structure["_refine_ls_number_reflns"]+"/"+ structure["_refine_ls_number_restraints"]+ "/"+ structure["_refine_ls_number_parameters"]  +" \\ \n")
    myfile.write("Goodness-of-fit                       &"+ structure["_refine_ls_restrained_S_all"] +" \\ \n")
    myfile.write("Final R indexes {[}I$geq$2$\\sigma$ (I){]} & R$_1$ = "+ structure["_refine_ls_R_factor_gt"]+ ", wR$_2$ = "+ structure["_refine_ls_wR_factor_gt"] + "\\\ \n")
    myfile.write("Final R indexes {[}all data{]}               & R$_1$ = "+ structure["_refine_ls_R_factor_all"]+ ", wR$_2$ = "+ structure["_refine_ls_wR_factor_ref"] + " \\\ \n")
    myfile.write("Largest diff. peak/hole / e \AA$^-3$            & "+ structure["_refine_diff_density_max"] + "/ "+ structure["_refine_diff_density_min"] +" \\* \\bottomrule \n")
    myfile.write("\\end{longtable}"

def atom_tables(filename, structure):
  atom_label = structure["_atom_site_label"]
  atom_x = structure["_atom_site_fract_x"]
  atom_y = structure["_atom_site_fract_y"]
  atom_z = structure["_atom_site_fract_z"]
  atom_U = structure["_atom_site_U_iso_or_equiv"]
  data_iso = pd.DataFrame({"Atom": atom_label, "x": atom_x , "y": atom_y, "z": atom_z, "U": atom_U})

  with open(filename+"_atom_coordinates_iso.tex", "w") as myfile:
    myfile.write(data_iso.to_latex(index=False))

  atom_anis = structure["_atom_site_aniso_label"]
  atom_u11 = structure["_atom_site_aniso_U_11"]
  atom_u22 = structure["_atom_site_aniso_U_22"]
  atom_u33 = structure["_atom_site_aniso_U_33"]
  atom_u23 = structure["_atom_site_aniso_U_23"]
  atom_u13 = structure["_atom_site_aniso_U_13"]
  atom_u12 = structure["_atom_site_aniso_U_12"]
  data_anis = pd.DataFrame({"Atom": atom_anis, "U_{11}": atom_u11 , "U_{22}": atom_u22 , "U_{33}": atom_u33 , "U_{23}": atom_u23 , "U_{13}": atom_u13, "U_{12}": atom_u12})

  with open(filename+"_atom_anis.tex", "w") as myfile:
    myfile.write(data_anis.to_latex(index=False))
    
if __name__ == '__main__':
	ciftable_tex(sys.argv[2], structure(sys.argv[2]))
	atom_tables(sys.argv[2], structure(sys.argv[2]))

