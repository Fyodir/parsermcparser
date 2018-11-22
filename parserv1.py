# An XML parser in Python for reading an LRG file in .xml format and exporting exons coordinates

# importing required modules
import xml.etree.ElementTree as ET
import sys

# read input filename from argument
fileName = sys.argv[1]  # type: xml_file

# Check file name is valid .xml
assert fileName[-4:] == '.xml', 'You have the wrong input file'

# Check no additional arguments provided on command line
assert len(sys.argv) < 3, "Too many arguments"

# Inputting XML file
tree = ET.parse(fileName)

root = tree.getroot()


# Function to return exon numbers from lrg file into a dictionary
def exon_counter(root):
    exon_counter_dict = {}
    for child in root[0]:
        if child.tag == "transcript" and child.attrib['name'] == "t1":  # only takes the first transcirpt
            for i in child:
                if i.tag == "exon":
                    if (str(i.attrib))[-4] == "'":
                        exon_counter_dict["Exon: " + (str(i.attrib))[-3:-2]] = None
                    elif (str(i.attrib))[-5] == "'":
                        exon_counter_dict["Exon: " + (str(i.attrib))[-4:-2]] = None
                    elif (str(i.attrib))[-6] == "'":
                        exon_counter_dict["Exon: " + (str(i.attrib))[-5:-2]] = None
                    else:
                        print("Something is wrong with the XML file")
    return(exon_counter_dict)


print(exon_counter(root))
