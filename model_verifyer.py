"""
    *******************************************************************************************************************************
    *                                                                                                                             *
    *  This module reads excel files containing reaction stoichiometries (eahc excel file contains one reaction stoichiometry)    *
    *  And reaction copounds and their chemical compositions for each reaction in separate excel files                            *
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



path1 = '/Users/makb047/Library/CloudStorage/OneDrive-TheUniversityofAuckland/UoA/Codes/Stoichiometric-Matrix/docs/Reaction1.xlsx'

r1 = xls.read_excel_reactions(path1)

path2 = '/Users/makb047/Library/CloudStorage/OneDrive-TheUniversityofAuckland/UoA/Codes/Stoichiometric-Matrix/docs/Reaction2.xlsx'

r2 = xls.read_excel_reactions(path2)

path3 = '/Users/makb047/Library/CloudStorage/OneDrive-TheUniversityofAuckland/UoA/Codes/Stoichiometric-Matrix/docs/Reaction3.xlsx'

r3 = xls.read_excel_reactions(path3)

reactions = [r1, r2, r3]

com_indices, sto_mat = smb.sto_mat_builder( reactions )



r1_coms_path = '/Users/makb047/Library/CloudStorage/OneDrive-TheUniversityofAuckland/UoA/Codes/Stoichiometric-Matrix/docs/Reaction1-Compounds.xlsx'

r2_coms_path = '/Users/makb047/Library/CloudStorage/OneDrive-TheUniversityofAuckland/UoA/Codes/Stoichiometric-Matrix/docs/Reaction2-Compounds.xlsx'

r3_coms_path = '/Users/makb047/Library/CloudStorage/OneDrive-TheUniversityofAuckland/UoA/Codes/Stoichiometric-Matrix/docs/Reaction3-Compounds.xlsx'

r1_comp = xls.read_excel(r1_coms_path)

r2_comp = xls.read_excel(r2_coms_path)

r3_comp = xls.read_excel(r3_coms_path)

reactions_comp = [ r1_comp, r2_comp, r3_comp ]

compounds, ele_indices = comex.comp_extract(reactions_comp)

element_matrix = emb.elemental_matrix_builder( com_indices, ele_indices, compounds )

sym_list = list(com_indices.keys())

rate_mat = rmb.rate_mat_builder (sym_list)



# Calling the function
vf.verification(sto_mat, element_matrix, rate_mat )