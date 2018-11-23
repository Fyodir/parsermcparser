# An XML parser in Python for reading an LRG file in .xml format and exporting exons coordinates

# importing required modules
import xml.etree.ElementTree as ET
import sys
import numpy as np
import os
import re

# read input filename from argument
fileName = sys.argv[1]  # type: xml_file

# Check file name is valid .xml
assert fileName[-4:] == '.xml', 'You have the wrong input file'

# Check no additional arguments provided on command line
assert len(sys.argv) < 3, "Too many arguments"

# Inputting XML file
tree = ET.parse(fileName)

root = tree.getroot()

for gene_name in tree.findall('.//lrg_locus'):
    gene = gene_name.text

# function to return exon numbers of lrg file in a list
def exon_num(root):
    exon_num_list = []
    for child in root[0]:
        if child.tag == "transcript" and child.attrib['name'] == "t1":  # only takes the first transcirpt
            for i in child:
                if i.tag == "exon":
                    if (str(i.attrib))[-4] == "'":
                        exon_num_list.append((str(i.attrib))[-3:-2]) # returns exons 1 - 9
                    elif (str(i.attrib))[-5] == "'":
                        exon_num_list.append((str(i.attrib))[-4:-2]) # returns exons 10 - 99
                    elif (str(i.attrib))[-6] == "'":
                        exon_num_list.append((str(i.attrib))[-5:-2]) # returns exons 100 - 999
                    else:
                        print("Something is wrong with the XML file")
    return(exon_num_list)


exon_num_var = np.asarray(exon_num(root))


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
start_list = np.asarray(start_list[0])
end_list = np.asarray(end_list[0]) # Pulls first entry from nested lists

header = "Exon\tStart\tEnd\n" # headers for output text file


# writing output file named by gene name, including exon number & LRG coordinates & headers

with open('%s.bed' % gene, 'w+') as file_temp:
    file_temp.write(header)

with open('%s.bed' % gene, 'a') as file_temp:
    for (exon_num_var, start_list, end_list) in zip(exon_num_var, start_list, end_list):
        file_temp.write("{0}\t\t{1}\t{2}\n".format(exon_num_var, start_list, end_list))




#print(gene)
#print(exon_num_var)
#print(start_list)
#print(end_list)


#a = [exon_num_var, start_list, end_list]
#nparray = (np.asarray(a))
#np.savetxt("test.txt", np.c_[exon_num_var, start_list, end_list], delimiter='\t')
