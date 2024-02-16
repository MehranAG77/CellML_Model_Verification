
"""
    *******************************************************************************************************************************
    *                                                                                                                             *
    *  This module reads excel files containing reaction stoichiometries (each excel file contains one reaction stoichiometry)    *
    *  And reaction copounds from a CellML file                                                                                   *
    *  Then using the modules provided constructs stoichiometric and elemental matrices and verifies the model                    *
    *                                                                                                                             *
    *******************************************************************************************************************************
"""

# Importing external or built-in packages
import numpy as np
import sympy as sp

# Importing internal packages
import excel_read as xls
import compounds_extractor as comex
import stochio_matrix_builer as smb
import elemental_matrix_builder as emb
import rate_matrix_builder as rmb
import verification as vf
import CellML_reader as cmlr
import compound_element_sorter as ces
import equation_builder as eb


cellml_file_dir = './docs/reactions_set.cellml'
cellml_file = './docs/reactions_set.cellml'
cellml_strict_mode = False

component = cmlr.CellML_reader( cellml_file, cellml_file_dir, cellml_strict_mode )

eb.equation_builder( component )

ele_indices, compound_indices, sym_list, compounds, sto_mat = ces.cellml_compound_element_sorter ( component )

element_matrix = emb.elemental_matrix_builder( compound_indices, ele_indices, compounds )

rate_mat = rmb.rate_mat_builder ( sym_list )

# Calling the function
vf.verification( sto_mat, element_matrix, rate_mat )
