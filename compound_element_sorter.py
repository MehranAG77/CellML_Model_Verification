
"""
CellML reader now gives a dictionary containing variable names as keys and variable ids as values. Variable id is like 'id_chebicode', e.g., id_18654. so we nned to extract the chebicode
and give to chemparse to parse it into its constituents and return it as a dictionary.
Then we put these chemicals into compounds dictionary where keys are chemical formula and values are dictionaries containing its constituents and their corresponding quantity
Elements dictionary which contains indices of elements is also constructed and returned
In the latest update (April 16, 2024), variable annotation can include the composition of the compound which is given by its name to CellML: 
'va_mekppmapk-mek.mapk': This annotation containes two parts separated by an underline(_). 2nd part contains the compound name and its composition. The minus sign(-) separates component name and composition.
So the first part (mekppmapk) is the component name. The 2nd part (mek.mapk) shows that it is composed of 'mek' and 'mapk'.
"""

import numpy as np
import chemparse as chp

from libcellml import Parser, Printer, Validator, cellmlElementTypeAsString
# Importing internal packages
from utilities import print_model
import cellml

import chebi_fetch as chf

import sys

# This function determines if the characters used in a variable are all digits or alphabets
def all_digits( word ):

    return all( char.isdigit() for char in word )


# This function sorts the variables of the CellMl file into the corresponding lists
def variable_sorter( components, all_vars = 'Y' ):

    variables = []

    coefficients = []

    rate_constants = []

    rates = []

    for component in components:

        number_of_variables = component.variableCount()

        if ( all_vars == "Y" or all_vars == "y" ): 
            
            for v in range( 0, number_of_variables ):
                
                if component.variable(v).id():      # There are some variables in CellML like time (t) that I don't want to import to python,so I left it without id and here I can import those variables with id that I need here

                    id = component.variable(v).id()

                    identifier = id.split('_')[0]

                    if identifier == 'va': variables.append( component.variable(v) )    # Since we have two different types of parameters in CellML, I put them in different lists
                    elif identifier == 'co': coefficients.append( component.variable(v) )
                    elif identifier == 'rc': rate_constants.append( component.variable(v) )
                    elif identifier == 'ra': rates.append( component.variable(v) )

        else:
                    
            for v in range( 0, number_of_variables ):
                
                if component.variable(v).id():      # There are some variables in CellML like time (t) that I don't want to import to python,so I left it without id and here I can import those variables with id that I need here

                    id = component.variable(v).id()

                    identifier = id.split('_')[0]

                    if identifier == 'va': variables.append( component.variable(v) ) 
                    elif identifier == 'co': coefficients.append( component.variable(v) )
                    

    if ( all_vars == "Y" or all_vars == "y" ): 

        return variables, coefficients, rates, rate_constants
    
    else:

        return variables, coefficients

def variable_name_mapper( components ):

    '''
    This function receives the components of a CellML file and maps the chebi codes to the variable name,
    and returns a dictionary mapping chebi code to the CellML variable name
    chebi_toCellML = { 'C2H6O12': x_glu, ... }
    '''

    chebi_to_CellML = {}

    chebi_initialvalues = {}

    variables, _ = variable_sorter( components , 'n')

    for variable in variables:      # After putting the parameters in their corresponding list, I get their ChEBI code and build the matrices
            
        name = variable.name()

        intial_value = variable.initialValue()

        chebi_code =  variable.id().split('_')[1]

        if all_digits( chebi_code ):

            compound, _ = chf.chebi_comp_parser( chebi_code )

        else:

            compound = chebi_code.split('-')[0]

        chebi_initialvalues[compound] = intial_value

        chebi_to_CellML[compound] = name

    return chebi_to_CellML, chebi_initialvalues



def cellml_compound_element_sorter ( components ):

    compounds = {}  # This dictionary will store the compounds and their chemical composition

    reaction_indices = {} # This dictionary keeps an index for each reaction that is used in constructing stoichiometric matrix ==> {'r1': 0, 'r2': 1, ...} "It can be unordered numbers for reactions or any alphabets ==> {'r10': 0, ...}"

    compound_indices = {} # This dictionary keeps an index for each compound to be used as a reference to construct stoichiometric matrix ==> {'C6H12O6': 0, 'CO2': 1, ...} so we know that row 0 in matrix is for C6H12O6

    element_indices ={}    # This dictionary contains elements as keys and assigns a number for each in order to construct the elemental matrix as values ==> {'C': 0, 'H': 1, 'O': 2}

    compound_index = 0     # This variable keeps the ordered number for assigning any compound a number

    element_index = 0     # This index keeps track of the order each new element will be put in the elemental matrix

    reaction_index = 0     # This index keeps track of the order for reactions to be used as a reference for the column in stoichiometric matrix

    variables = []      # This list contains the variables of cellml file [In CellML file, I defined two types of parameters, first type is variable which can change in the model and the second type is coefficient which is constant]

    coefficients = []       # The difference between variables and coefficients is identified by the first part of their id which is like "va_12345_r1", first part indicates that this parameter is a variable

    symbols_list = []       # This list contains 

    number_of_components = len( components )

    variables, coefficients = variable_sorter( components , 'n')

    for variable in variables:      # After putting the parameters in their corresponding list, I get their ChEBI code and build the matrices
            
        name = variable.name()

        chebi_code =  variable.id().split('_')[1]

        if all_digits( chebi_code ):

            compound, composition = chf.chebi_comp_parser( chebi_code )

        else:
            
            if ( len(chebi_code.split('-')) == 1 ):
                compound = chebi_code
                composition={ chebi_code: 1 }
            
            else:

                compound = chebi_code.split('-')[0]
                formula = chebi_code.split('-')[1].split('.')

                composition={}
                for element in formula:
                    composition[element] = 1

        if compound not in compounds:
                
            compounds[compound] = composition

        if compound not in compound_indices:
                
            compound_indices[compound] = compound_index
            compound_index += 1
            
        if name not in symbols_list:

            symbols_list.append( name )

        for element in composition:
                        
            if element not in element_indices:

                element_indices[element] = element_index
                element_index += 1

    for coefficient in coefficients:
            
        chebi_code = coefficient.id().split('_')[1]

        reaction_number = coefficient.id().split('_')[2]

        if reaction_number not in reaction_indices:
            reaction_indices[reaction_number] = reaction_index
            reaction_index += 1

    columns = len( reaction_indices ) # Number of columns will be equal to number of reactions in the Stoichiometric matrix

    rows = len( compound_indices ) # Number of rows will be equal to number of compounds in the Stoichiometric matrix

    # print("The reactions are: {l1}".format(l1=list(reaction_indices.keys())))
    # print("The Compounds are: {l2}".format(l2=list(compound_indices.keys())))

    stoichiometric_matrix = np.zeros(( rows, columns), dtype = int)

    for coefficient in coefficients:

        chebi_code = coefficient.id().split('_')[1]

        reaction_number = coefficient.id().split('_')[2]

        if all_digits( chebi_code ):

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
        

    return element_indices, compound_indices, symbols_list, compounds, stoichiometric_matrix, reaction_indices