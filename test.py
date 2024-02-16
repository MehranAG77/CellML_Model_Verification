from lxml import etree
import xml.dom.minidom
from xml_flatten import xml_flatten
import re
import xml.etree.ElementTree as ET

def normalize_mathml(mathml_string):
    root = ET.fromstring(mathml_string)
    normalized_string = ET.tostring(root, encoding='unicode', method = 'xml', short_empty_elements=False )
    return normalized_string

text2 = """
<apply xmlns:mml=\"http://www.w3.org/1998/Math/MathML\"><eq/><apply><diff/><bvar><ci>t</ci></bvar><apply><x_ethanol/><ci>t</ci></apply></apply><apply><times/><ci><mml:msub><mml:mi>ethanol</mml:mi><mml:mi>r3</mml:mi></mml:msub></ci><ci><mml:msub><mml:mi>v</mml:mi><mml:mi>3</mml:mi></mml:msub></ci></apply></apply>
"""

text1 = """
<apply>
                <eq/>
                <apply>
                    <diff/>
                    <bvar>
                        <ci>t</ci>
                    </bvar>
                    <ci>x_ethanol</ci>
                </apply>
                <apply>
                    <times/>
                    <ci>ethanol_r3</ci>
                    <ci>v_3</ci>
                </apply>
            </apply>
"""

def prettify_xml(one_row_mathml):
    dom = xml.dom.minidom.parseString(one_row_mathml)

    # Pretty-print the XML with indentation
    tidy_xml = dom.toprettyxml(indent="    ")

    return(tidy_xml)


def clean_xml_indent(xml_string):
    # Use regular expression to remove leading whitespace and indentation
    cleaned_xml = re.sub(r'\s+', ' ', xml_string)
    cleaned_xml = re.sub(r'>\s+<', '><', cleaned_xml)

    return cleaned_xml

def remove_first_row(input_string):
    # Split the string into lines
    lines = input_string.splitlines()

    # Check if there is more than one line
    if len(lines) > 1:
        # Remove the first line
        lines.pop(0)

    # Join the lines back into a string
    result_string = '\n'.join(lines)

    return result_string

text2 = prettify_xml(text2)

text2 = remove_first_row( text2 )

#print(text1)

root1 = etree.fromstring(text1)
root2 = etree.fromstring(text2)


mylist= root2.xpath('//ci')


expression_to_find = []

replacement_expression = []

segment_to_replace = "<ci xmlns:mml=\"http://www.w3.org/1998/Math/MathML\">"
replacement_string = "<ci>"

for item in mylist:
    if len(item) > 0:
        mystring = etree.tostring(item, pretty_print=True, encoding=str)
        #cleaned_xml = clean_xml_indent(mystring)
        #prettified_cleaned_xml = prettify_xml(cleaned_xml)
        #my_new_string = prettified_cleaned_xml.replace(segment_to_replace, replacement_string)
        #my_newest_string = remove_first_row(my_new_string)
        expression_to_find.append(mystring)
        alist=[]
        for item1 in item.iter():
            variable = item1.text.lstrip('\n')
            if variable.isspace():
                continue
            alist.append(variable)
        new_variable = '<ci>'+alist[0]+'_'+alist[1]+'</ci>'
        replacement_expression.append(new_variable)


def replacement( mathml_script, expression_to_find, replacement_expression ):

    # Parse the MathML script
    root = etree.fromstring(mathml_script)

    # Find and replace the expression
    for elem in root.iter():
        # Check if the current element matches the expression to find
        #print('\nElement is\n', elem)
        #print('\nString Element is\n', etree.tostring(elem, encoding="unicode", method='xml', xml_declaration=None, pretty_print=True))
        #print('\nExpression to find is\n', expression_to_find.strip())
        if normalize_mathml(etree.tostring(elem, encoding="unicode")) == normalize_mathml(expression_to_find.strip()):
            # Replace the element with the replacement expression
            parent = elem.getparent()
            index = parent.index(elem)
            parent.remove(elem)

            replacement_elem = etree.fromstring(replacement_expression)
            parent.insert(index, replacement_elem)
            parent[index].tail = "\n" + "\t"

    # Print the modified MathML script
    modified_mathml_script = etree.tostring(root, encoding="unicode", pretty_print=True)
    print("Modified MathML Script:")
    print( modified_mathml_script )
    return modified_mathml_script


for i in range(0, len(expression_to_find)):

   text2 = replacement( text2, expression_to_find[i], replacement_expression[i] )