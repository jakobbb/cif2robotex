from pymatgen.io.cif import CifParser
from pymatgen import Structur

def cif_import(filename):
  parser = CifParser(filename + ".cif")
  structure = parser.as_dict()
  structure = structure[list(structure.keys())[0]]
  return structure

def cif_import(filename):
  parser = CifParser(filename + ".cif")
  structure = parser.as_dict()
  structure = structure[list(structure.keys())[0]]
  return structure

import pandas as pd
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

