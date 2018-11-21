# An XML parser in Python for reading an LRG file in .xml format and exporting exons coordinates

# importing required modules
import xml.etree.ElementTree as ET
import sys

# read input filename from argument
fileName = sys.argv[1]  # type: xml file

#Check file name is valid .xml
assert fileName[-4:] == '.xml', 'You have the wrong input file'

# Inputting XML file
tree = ET.parse(fileName)

root = tree.getroot()


# prints exon numbers of lrg file
for child in root[0][8]:
    if child.tag == "exon":
        if (str(child.attrib))[-4] == "'":
            print ("Exon: " + (str(child.attrib))[-3:-2])
        else:
            print("Exon: " + (str(child.attrib))[-4:-2])
