#Parsing module
import xml.etree.ElementTree as ET

#Inputting XML file
tree = ET.parse('inputlrg')

root = tree.getroot()
