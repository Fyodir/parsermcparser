# An XML parser in Python for reading an LRG file in .xml format and exporting exons coordinates

# importing required modules
import xml.etree.ElementTree as ET
import sys
import numpy as np

# read input filename from argument
fileName = sys.argv[1]  # type: xml_file

# Check file name is valid .xml
assert fileName[-4:] == '.xml', 'You have the wrong input file'

# Check no additional arguments provided on command line
assert len(sys.argv) < 3, "Too many arguments"

# Inputting XML file
tree = ET.parse(fileName)

root = tree.getroot()


# function to return exon numbers of lrg file in a list
def exon_num(root):
    exon_num_list = []
    for child in root[0]:
        if child.tag == "transcript" and child.attrib['name'] == "t1":  # only takes the first transcirpt
            for i in child:
                if i.tag == "exon":
                    if (str(i.attrib))[-4] == "'":
                        exon_num_list.append("Exon: " + (str(i.attrib))[-3:-2]) # returns exons 1 - 9
                    elif (str(i.attrib))[-5] == "'":
                        exon_num_list.append("Exon: " + (str(i.attrib))[-4:-2]) # returns exons 10 - 99
                    elif (str(i.attrib))[-6] == "'":
                        exon_num_list.append("Exon: " + (str(i.attrib))[-5:-2]) # returns exons 100 - 999
                    else:
                        print("Something is wrong with the XML file")
    return(exon_num_list)


exon_num_var = exon_num(root)


# function to return two lists for start and end coordinates of exons respectively
def exon_coord(root):
    exon_start_list = []
    exon_end_list = []
    for child in root[0]:
        if child.tag == "transcript" and child.attrib['name'] == "t1":  # only takes the first transcirpt
            for exon in child:
                if exon.tag == "exon":
                   exon_start_list.append(exon[0].attrib["start"])
            for exon in child:
                if exon.tag == "exon":
                   exon_end_list.append(exon[0].attrib["end"])
    return exon_start_list, exon_end_list

start_list, end_list = map(list, zip(exon_coord(root)))
start_list = start_list[0]
end_list = end_list[0]



#print(exon_num_var)
#print(start_list)
#print(end_list)
#print(exon_num_var)


array = [exon_num_var, start_list, end_list]
print(np.asarray(array))


