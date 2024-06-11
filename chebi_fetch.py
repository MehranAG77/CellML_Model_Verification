# Importing external packages
import libchebipy as chb
import chemparse as chp


# chebi_code is a string which contains ChEBI code like '46683'
def chebi_comp_parser( chebi_code ):

    """
    This function receives a string which includes ChEBI code for the compound and uses EBI API to search for the compounds chemical composition and fetches it
    Then using the chemparse package decomposes the chemical formula to its elements and returns it as a dictionary
    It should be mentioned that some ChEBI compounds don't have chemical formula registered for them and some have two different chemical compositions assigned to them
    """

    chebi_entity = chb.ChebiEntity(chebi_code)
    name = chebi_entity.get_name()

    parsed_compound={}

    if __name__ == '__main__':
        print(f"\nThe name of this compound is: {name}")

    formulae = chebi_entity.get_formulae() # There are two methods to fetch chemical composition of a compound: get_formulae is used when more than one composition is registerd for the compound
                                            # It returns a list containing 'Formula' classes and Getting list's length will show us that if there is one or more chemical compositions registered
    if len(formulae) == 0: # Sometimes no chemical formula is registered so the length of the list will be zero

        if __name__ == '__main__':
            print("No formula is registered for this compound")

        parsed_compound = None
        formula = None

    elif len(formulae) == 1:

        formula = chebi_entity.get_formula()
        parsed_compound = chp.parse_formula(formula)
        if __name__ == '__main__':
            print( "Chemical formula for this copmound is: {v}".format( v = formula ) )
            print(chp.parse_formula(formula))

    else:

        formula = formulae[1].get_formula()
        parsed_compound = chp.parse_formula(formula)

        if __name__ == '__main__':
            print( "Chemical formula for this copmound is: {v}".format( v = formula ) )
            print(chp.parse_formula(formula))

    return formula, parsed_compound



def chebi_formula ( chebi_code ):
    
    """
    This function receives a string which includes ChEBI code for the compound and uses EBI API to search for the compound's chemical composition and fetches it
    It should be mentioned that some ChEBI compounds don't have chemical formula registered for them and some have two different chemical compositions assigned to them
    """

    chebi_entity = chb.ChebiEntity(chebi_code)

    formulae = chebi_entity.get_formulae() # There are two methods to fetch chemical composition of a compound: get_formulae is used when more than one composition is registerd for the compound
                                            # It returns a list containing 'Formula' classes and Getting list's length will show us that if there is one or more chemical compositions registered
    if len(formulae) == 0: # Sometimes no chemical formula is registered so the length of the list will be zero

        if __name__ == '__main__':
            print("No formula is registered for this compound")

        formula = None

    elif len(formulae) == 1:

        formula = chebi_entity.get_formula()

        if __name__ == '__main__':
            print( "Chemical formula for this copmound is: {v}".format( v = formula ) )

    else:

        formula = formulae[1].get_formula()

        if __name__ == '__main__':
            print( "Chemical formula for this copmound is: {v}".format( v = formula ) )

    return formula






if __name__ == '__main__':

    chebi_code = '18367'

    formula, parsed = chebi_comp_parser( chebi_code )

    print(formula, parsed)