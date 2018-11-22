# An XML parser in Python for reading an LRG file in .xml format and exporting exons coordinates

# importing required modules
import xml.etree.ElementTree as ET
import sys

# read input filename from argument
fileName = sys.argv[1]  # type: xml_file

#Check file name is valid .xml
assert fileName[-4:] == '.xml', 'You have the wrong input file'

# Inputting XML file
tree = ET.parse(fileName)

root = tree.getroot()


# prints exon numbers of lrg file
for child in root[0]:
    if child.tag == "transcript":
        for i in child:
            if i.tag == "exon":
                if (str(i.attrib))[-4] == "'":
                    print ("Exon: " + (str(i.attrib))[-3:-2])
                elif (str(i.attrib))[-5] == "'":
                    print("Exon: " + (str(i.attrib))[-4:-2])
                elif (str(i.attrib))[-6] == "'":
                    print ("Exon: " + (str(i.attrib))[-5:-2])

