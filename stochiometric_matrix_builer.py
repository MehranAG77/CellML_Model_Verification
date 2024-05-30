"""
This function receives a list containing a numpy array of compounds and their corresponding stoichiometric coefficients
And returns a dictionary containing compounds and their indices showing which row they are assigned in the stoichiometric matrix
And also returns the stoichiometic matrix (numpy array)
"""

# Importing Built-in packages
import sys


# Importing external or built-in packages
import numpy as np


# Importing internal packages
import compound_element_sorter as ces
import chebi_fetch as chf




def stoichiometric_matrix_builder( reaction_indices, compound_indices, coefficients, bc_coefficients ):

    columns = len( reaction_indices ) # Number of columns will be equal to number of reactions in the Stoichiometric matrix

    rows = len( compound_indices ) # Number of rows will be equal to number of compounds in the Stoichiometric matrix

    stoichiometric_matrix = np.zeros(( rows, columns), dtype = int)

    for coefficient in coefficients:

        chebi_code = coefficient.id().split('_')[1]

        reaction_number = coefficient.id().split('_')[2]

        if ces.all_digits( chebi_code ):

            formula = chf.chebi_formula( chebi_code )

        else:

            formula = chebi_code

        try:

            row = compound_indices[formula]

        except KeyError as KE:

            print("The stoichiometric coefficient for compound {v} cannot be found in your coefficients of CellML".format( v = KE ))
            sys.exit("Exiting due to an error\nModify CellML file and check to see if you have this compound or errors in a ChEBI code in reaction {r}".format( r = reaction_number))

        column = reaction_indices[reaction_number]

        stoichiometric_matrix[ row, column ] = coefficient.initialValue()



    # Constructing the rows and columns for the boundary conditions in the stoichiometric matrix
    for bc_coefficient in bc_coefficients.keys():

        compound = bc_coefficient.split('-')[0]

        reaction_number = bc_coefficient.split('-')[1]

        try:

            row = compound_indices[compound]

        except KeyError as KE:

            print("The stoichiometric coefficient for compound {v} cannot be found in your coefficients of CellML".format( v = KE ))
            sys.exit("Exiting due to an error\nModify CellML file and check to see if you have this compound or errors in a ChEBI code in reaction {r}".format( r = reaction_number))

        column = reaction_indices[reaction_number]

        stoichiometric_matrix[ row, column ] = bc_coefficients[bc_coefficient]
        

    return stoichiometric_matrix