

"""
This function takes a list containing the dataframes (as a list variable) for reaction compound elements
"""


# Importing external and built-in packages
import numpy as np
import pandas as pd

# Importing internal packages
import excel_read as xls


def comp_extract( coms ):

    compounds = {}   # This dictionary will contain compounds as keys and a dictionary (containing compound elements and their respective subscript) as values

    ele_dict ={}    # This dictionary contains elements as keys and assigns a number for each in order to construct the elemental matrix as values

    reactions_no = len(coms)    # At first, the number of reactions are fetched by knowing the length of the list

    m_index = 0 # This index will keep track of the number of reaction being explored

    while m_index < reactions_no:

        counter = 0 # This value keeps track of the number of compounds in the reaction

        e_index = 0 # This index keeps track of the order each new element will be put in the elemental matrix

        while counter < len(coms[m_index].columns):


            com = coms[m_index].columns[counter][0] # The name of compound is stored in com

            counter = counter + 2   # Since each compound has two columns, counter is increased by 2 to reach the other compound in the reaction in the next iteration

            if com in compounds:    # If the compound exists in the dictionary so it is decoded and no need to do it again
                pass
            else:   # If the compound is new, we need to look into it

                index = 0   # This index keeps track of the element number
                mydict = {} # A dictionary is created to store elements

                while index < len(coms[m_index][com]['E']): #Sweeping rows for elements of the compound
            
                    if pd.isnull(coms[m_index][com]['E'][index]):
                        pass
                    else:
                        mydict.update({coms[m_index][com]['E'][index]:coms[m_index][com]['C'][index]})
                    index += 1
                compounds[com] = mydict


                for ele in compounds[com]:
                    if ele in ele_dict:
                        pass
                    else:
                        ele_dict[ele] = e_index
                        e_index += 1

        m_index += 1

    return compounds, ele_dict