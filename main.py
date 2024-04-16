
"""
    *******************************************************************************************************************************
    *                                                                                                                             *
    *  This module reads CellML files containing reaction stoichiometries and equations                                           *
    *                                                                                                                             *
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
import stochiometric_matrix_builer as smb
import elemental_matrix_builder as emb
import rate_matrix_builder as rmb
import verification as vf
import CellML_reader as cmlr
import compound_element_sorter as ces
import equation_builder as eb


cellml_file_dir = './docs/aguda_b_1999.cellml'
cellml_file = './docs/aguda_b_1999.cellml'
cellml_strict_mode = False

components = cmlr.CellML_reader( cellml_file, cellml_file_dir, cellml_strict_mode )

# eb.equation_builder( components )

element_indices, compound_indices, symbols_list, compounds, stoichiometric_matrix = ces.cellml_compound_element_sorter ( components )

element_matrix = emb.elemental_matrix_builder( compound_indices, element_indices, compounds )

rate_matrix = rmb.rate_matrix_builder ( symbols_list )

# Calling the function
vf.verification( stoichiometric_matrix, element_matrix, rate_matrix )
