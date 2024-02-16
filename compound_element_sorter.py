
"""
CellML reader now gives a dictionary containing variable names as keys and variable ids as values. Variable id is like 'id_chebicode', e.g., id_18654. so we nned to extract the chebicode
and give to chemparse to parse it into its constituents and return it as a dictionary.
Then we put these chemicals into compounds dictionary where keys are chemical formula and values are dictionaries containing its constituents and their corresponding quantity
Elements dictionary which contains indices of elements is also constructed and returned
"""

import numpy as np
import chemparse as chp

from libcellml import Parser, Printer, Validator, cellmlElementTypeAsString
# Importing internal packages
from utilities import print_model
import cellml

import chebi_fetch as chf

import sys


def cellml_compound_element_sorter ( component ):

    compounds = {}  # This dictionary will strore the compounds and their chemical composition

    reaction_indices = {} # This dictionary keeps an index for each reaction that is used in constructing stoichiometric matrix ==> {'r1': 0, 'r2': 1, ...} "It can be unordered numbers for reactions or any alphabets ==> {'r10': 0, ...}"

    compound_indices = {} # This dictionary keeps an index for each compound to be used as a reference to construct stoichiometric matrix ==> {'C6H12O6': 0, 'CO2': 1, ...} so we know that row 0 in matrix is for C6H12O6

    ele_indices ={}    # This dictionary contains elements as keys and assigns a number for each in order to construct the elemental matrix as values ==> {'C': 0, 'H': 1, 'O': 2}

    c_index = 0     # This variable keeps the ordered number for assigning any compound a number

    e_index = 0     # This index keeps track of the order each new element will be put in the elemental matrix

    r_index = 0     # This index keeps track of the order for reactions to be used as a reference for the column in stoichiometric matrix

    number_of_variables = component.variableCount()

    variables = []      # This list contains the variables of cellml file [In CellML file, I defined two types of parameters, first type is variable which can change in the model and the second type is coefficient which is constant]

    coefficients = []       # The difference between variables and coefficients is identified by the first part of their id which is like "va_12345_r1", first part indicates that this parameter is a variable

    sym_list = []       # This list contains 

    for v in range( 0, number_of_variables ):

        if component.variable(v).id():      # There are some variables in CellML like time (t) that I don't want to import to python,so I left it without id and here I can import those variables with id that I need here

            id = component.variable(v).id()

            identifier = id.split('_')[0]

            if identifier == 'va': variables.append( component.variable(v) )    # Since we have two different types of parameters in CellML, I put them in different lists
            elif identifier == 'co': coefficients.append( component.variable(v) )

    for variable in variables:      # After putting the parameters in their corresponding list, I get their ChEBI code and build the matrices
        
        name = variable.name()

        chebi_code =  variable.id().split('_')[1]

        formula, composition = chf.chebi_comp_parser( chebi_code )

        if formula not in compounds:
             
            compounds[formula] = composition

        if formula not in compound_indices:
             
             compound_indices[formula] = c_index
             c_index += 1
        
        if name not in sym_list:

            sym_list.append( name )

        for ele in composition:
                    
            if ele not in ele_indices:

                ele_indices[ele] = e_index
                e_index += 1

    for coefficient in coefficients:
         
        chebi_code = coefficient.id().split('_')[1]

        re_number = coefficient.id().split('_')[2]

        if re_number not in reaction_indices:
             reaction_indices[re_number] = r_index
             r_index += 1

    columns = len( reaction_indices ) # Number of columns will be equal to number of reactions in the Stoichiometric matrix

    rows = len( compound_indices ) # Number of rows will be equal to number of compounds in the Stoichiometric matrix

    sto_mat = np.zeros(( rows, columns), dtype = int)

    for coefficient in coefficients:

        chebi_code = coefficient.id().split('_')[1]

        re_number = coefficient.id().split('_')[2]

        formula = chf.chebi_formula( chebi_code )

        try:

            row = compound_indices[formula]

        except KeyError as KE:

            print("The stoichiometric coefficient for compound {v} cannot be found in your coefficients of CellML".format( v = KE ))
            sys.exit("Exiting due to an error\nModify CellML file and check to see if you have this compound or errors in a ChEBI code in reaction {r}".format( r = re_number))

        column = reaction_indices[re_number]

        sto_mat[ row, column ] = coefficient.initialValue()
        

    return ele_indices, compound_indices, sym_list, compounds, sto_mat