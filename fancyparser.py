# An XML parser in Python for reading an LRG file in .xml format and exporting exons coordinates

# importing required modules
import xml.etree.ElementTree as ET
import sys
import numpy as np
from operator import sub
import time
import requests

input_lrg = input("Please enter LRG number: ")

#Assert to ensure only positive integers are entered by the user
assert input_lrg.isdigit(), "You must enter a positive integer"

url = 'http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_%s.xml' % input_lrg
r = requests.get(url, allow_redirects=True)
open('LRG_%s.xml' % input_lrg, 'wb').write(r.content)

fileName = open('LRG_%s.xml' % input_lrg, 'r')

# Inputting XML file
tree = ET.parse(fileName)
root = tree.getroot()

for gene_name in tree.findall('.//lrg_locus'):
    gene = 'LRG' + '_' + str(input_lrg) + "_" + gene_name.text

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

start_list_str, end_list_str = map(list, zip(exon_coord(root)))

# Converts start_list_str and end_list_str to integer values
lrg_start_list = []
lrg_end_list = []

for i in start_list_str[0]:
    lrg_start_list.append(int(i))

for i in end_list_str[0]:
    lrg_end_list.append(int(i))

# function to calculate exon lengths
def exon_len_func(a,b):
    exon_len = []
    exon_len_count = 0
    for i in a:
        exon_len.append(b[exon_len_count]-a[exon_len_count])
        exon_len_count += 1
    return exon_len

exon_len = exon_len_func(lrg_start_list, lrg_end_list)

for i in tree.findall('.//mapping'):
    if i.attrib["coord_system"] == "GRCh37.p13":
        chromosome = i.attrib["other_name"] # Pulls chromosome number from XML
        gene_chr_start = int(i[0].attrib["other_start"]) #for loop to obtain start coords of gene on GRCh37.p13
        gene_chr_end = int(i[0].attrib["other_end"]) #for loop to obtain end coords of gene on GRCh37.p13
        strand = int(i[0].attrib["strand"]) #identify between forward strand(1) and reverse strand(-1)

# Mapping LRG coords to chromosomal coordinates (FORWARD STRAND)
if strand == 1:
    chr_exon_start = []
    for coord in lrg_start_list:
        chr_exon_start.append(coord + gene_chr_start -1)

    chr_exon_end = []
    for coord in lrg_end_list:
        chr_exon_end.append(coord + gene_chr_start -1)

else: # Mapping of LRG coordinates to chromosomal locations (REVERSE STRAND)
    chr_exon_start = []
    for coord in lrg_start_list:
        chr_exon_start.append(gene_chr_end - coord + 1)
    chr_exon_end = []
    for coord in lrg_end_list:
        chr_exon_end.append(gene_chr_end - coord + 1)

if strand == 1:
    strand = "Forward Strand\n\n"
else:
    strand = "Reverse Strand\n\n"

# pulls chromosome number from input LRG_xml
chr_list = []
count = 0
while count < len(chr_exon_start):
    chr_list.append("chr" + chromosome)
    count += 1

# Creates a date/time stamp for creaton of BED file
date = time.strftime("File created: %d/%m/%Y  %H:%M:%S\n\n")

# Creation of output file named by gene name

# Includes date/time stamp, column headers, followed by various columns of data

header = "\tStart\t\tEnd\t\tExon\tLength\n" # headers for output text file

with open('%s.bed' % gene, 'w+') as file_temp:
    file_temp.write(date)
    file_temp.write(strand)
    file_temp.write(header)
    for (chr_list, chr_exon_start, chr_exon_end, exon_num_var, exon_len) in zip(chr_list, chr_exon_start, chr_exon_end, exon_num_var, exon_len):
        file_temp.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(chr_list, chr_exon_start, chr_exon_end, exon_num_var, exon_len))
