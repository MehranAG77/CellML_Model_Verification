

""" 
This function is used to build elemental matrix by having: 1- compound indices in stoichiometric matrix (by indices I mean the order in which compounds
are put in the columns of the stoichiometric matrix), 2- element indices as a dictionary 'ele_indices' so that the lements are put in order in the row and can be identified
by the function when assigning their number in the compound, 3- compounds and their chemical composition as a dictionary where maps keys of compounds as strings to
values as dictionaries containing keys as elements and values as their corresponding number in the chemical composition of the compound
and returns elemental matrix where columns represent compounds and rows represent elements


com_indices is like this: {'C2H6O':0, 'CO2':1}
ele_indices is like this: {'C':0, 'H':1, 'O':2}
compounds is like this: {'C2H6O':{'C':2, 'H':6, 'O':1}, 'CO2':{'C':1, 'O':2}}

elemental matrix:

         C2H6O    CO2
       _               _
    C |    1       1    |
    H |    6       0    |
    O |    1       2    |
       -               -

  """

# Importing external and built-in packages
import numpy as np
import sympy as sp

# Impoting internal or built-in packages
import sys


def elemental_matrix_builder ( compound_indices, element_indices, compounds ):
    
    """
    This function takes three dictionaris: 1- compound_indices: This dictionary contains compounds with their corresponding index.
    2- element_indices: This dictionary contains elements of the compounds in the model with their corresponding indices.
    3- compounds: This dictionary contains compounds with their compostion as another dictionary stored as a value in the current dictionary
    The function returns the Elemental matrix as anumpy array
    """

    rows = len( element_indices )   # The number of rows will be the number of elements we have so it is as the length of the element dictionary
    
    columns = len( compound_indices )    # The number of columns will be the number of compounds we have so it is as the length of the compounds dictionary

    elemental_matrix = np.zeros((rows, columns), dtype = int) # construction of an empty matrix containing only zeros so that we can assign a number for elements in the matrix and the rest will be zeros

    for compound, compound_index in compound_indices.items():

        column = compound_index
        try:
          for element, element_index in compounds[compound].items():

              row = element_indices[element]
              value = element_index
              elemental_matrix[row][column] = value
        except KeyError as KE:
            print("The compound {v} cannot be found in your variables of CellML".format( v = KE ))
            sys.exit("Exiting due to an error\nModify CellML file and add {v}".format( v = KE))

    return elemental_matrix