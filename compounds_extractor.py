
# Importing external and built-in packages
import numpy as np
import pandas as pd

# Importing internal packages
import excel_read as xls





def comp_extract( compounds_df ):

    """
    This function takes a list containing the dataframes (as a list variable) for reaction compound elements
    """

    compounds = {}   # This dictionary will contain compounds as keys and a dictionary (containing compound elements and their respective subscript) as values

    element_dictionary ={}    # This dictionary contains elements as keys and assigns a number for each in order to construct the elemental matrix as values

    reactions_number = len(compounds_df)    # At first, the number of reactions are fetched by knowing the length of the list

    m_index = 0 # This index will keep track of the number of reaction being explored

    while m_index < reactions_number:

        counter = 0 # This value keeps track of the number of compounds in the reaction

        element_index = 0 # This index keeps track of the order each new element will be put in the elemental matrix

        while counter < len(compounds_df[m_index].columns):


            compound = compounds_df[m_index].columns[counter][0] # The name of compound is stored in com

            counter = counter + 2   # Since each compound has two columns, counter is increased by 2 to reach the other compound in the reaction in the next iteration

            if compound in compounds:    # If the compound exists in the dictionary so it is decoded and no need to do it again
                pass
            else:   # If the compound is new, we need to look into it

                index = 0   # This index keeps track of the element number
                my_dictionary = {} # A dictionary is created to store elements

                while index < len(compounds_df[m_index][compound]['E']): #Sweeping rows for elements of the compound
            
                    if pd.isnull(compounds_df[m_index][compound]['E'][index]):
                        pass
                    else:
                        my_dictionary.update({compounds_df[m_index][compound]['E'][index]:compounds_df[m_index][compound]['C'][index]})
                    index += 1
                compounds[compound] = my_dictionary


                for element in compounds[compound]:
                    if element in element_dictionary:
                        pass
                    else:
                        element_dictionary[element] = element_index
                        element_index += 1

        m_index += 1

    return compounds, element_dictionary