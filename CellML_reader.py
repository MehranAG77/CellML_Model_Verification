
"""
This function receives three inputs: 1- cellml_file, which is the name of the cellml file as a string, 2- cellml_file_dir, which is the path as a string to the cellml file,
and 3- cellml_strict_mode, which is a boolean with True and False identifying CellML 2.0 and 1.1 versions, respectively.
The new version (updated on 12/12/2023) returns the class of a component (libcellml class) that can be used to find variables and coefficients by another module
"""


# Importing external packages
from libcellml import Parser, Printer, Validator, cellmlElementTypeAsString
import numpy as np
import chemparse as chp

# Importing internal packages
from utilities import print_model
import cellml
import chebi_fetch as chf
import sys



cellml_file_dir = '/Users/makb047/Library/CloudStorage/OneDrive-TheUniversityofAuckland/UoA/Codes/Stoichiometric-Matrix/docs/reactions_set.cellml'
cellml_file = '/Users/makb047/Library/CloudStorage/OneDrive-TheUniversityofAuckland/UoA/Codes/Stoichiometric-Matrix/docs/reactions_set.cellml'
cellml_strict_mode = False

def CellML_reader( cellml_file, cellml_file_dir, cellml_strict_mode ):        # strict_mode is True when CellML version is 2.0 and False when 1.0 and 1.1

    model = cellml.parse_model(cellml_file, cellml_strict_mode)
    if cellml.validate_model(model) > 0:
        exit(-1)

    importer = cellml.resolve_imports(model, cellml_file_dir, cellml_strict_mode)
    if model.hasUnresolvedImports():
        print("unresolved imports?")
        exit(-2)

    if cellml.validate_model(model) > 0:
        print('Validation issues found')
        exit(-3)

    print('Model was parsed, resolved, and validated without any issues.')

    # need a flattened model for analysing
    flat_model = cellml.flatten_model(model, importer)
    if cellml.validate_model(flat_model) > 0:
        print('Validation issues found in flattened model')
        exit(-4)

    print('Model was flattened without any issues.')

    component = model.component(0)
    component_name = component.name()
    component_id = component.id()

    return component

"""    compounds = {}  # This dictionary will strore the compounds and their chemical composition

    reaction_indices = {} # This dictionary keeps an index for each reaction that is used in constructing stoichiometric matrix

    compound_indices = {}

    ele_indices ={}    # This dictionary contains elements as keys and assigns a number for each in order to construct the elemental matrix as values

    c_index = 0

    e_index = 0 # This index keeps track of the order each new element will be put in the elemental matrix

    r_index = 0

    number_of_variables = component.variableCount()

    variables = []

    coefficients = []

    sym_list = []

    for v in range( 0, number_of_variables ):

        if component.variable(v).id():

            id = component.variable(v).id()

            identifier = id.split('_')[0]

            if identifier == 'va': variables.append( component.variable(v) )
            elif identifier == 'co': coefficients.append( component.variable(v) )

    for variable in variables:
        
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

    reactions_no = len( reaction_indices )

    rows = len( compound_indices )

    sto_mat = np.zeros(( rows, reactions_no), dtype = int)

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
        

    return ele_indices, compound_indices, sym_list, compounds, sto_mat """


if __name__ == '__main__':

    print(CellML_reader( cellml_file, cellml_file_dir, cellml_strict_mode ))