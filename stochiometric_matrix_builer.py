"""
This function receives a list containing a numpy array of compounds and their corresponding stoichiometric coefficients
And returns a dictionary containing compounds and their indices showing which row they are assigned in the stoichiometric matrix
And also returns the stoichiometic matrix (numpy array)
"""

# Importing external or built-in packages
import numpy as np

# Importing internal packages
import excel_read as xls


def stochiometric_matrix_builder( reactions ):

    compound_indices = {}

    counter = 0

    index = 0

    reactions_no = len(reactions)

    while index < len(reactions):

        for compound in reactions[index]:

            if compound[0] in compound_indices:
                pass
            else:
                compound_indices[compound[0]] = counter
                counter += 1
        index += 1

    sto_mat = np.zeros((counter, reactions_no), dtype = int)

    index = 0

    while index < len(reactions):

        for compound in reactions[index]:

            sto_mat[compound_indices[compound[0]], index]  = compound[1]
        
        index += 1

    return compound_indices, sto_mat

#print(sto_mat)