'''
An XML parser in Python for reading an LRG file in .xml format and exporting exons coordinates
Authors: Andrew Smith and Jethro Rainford
Development start date: 27th November 2018
Usage: refer to the README file
'''

# Importing required modules
import xml.etree.ElementTree as ET
import time
import requests


# Function to use user input to pull LRG XML from  LRG website and create local XML file for parsing
def lrg_input(input_lrg):
    # assert to ensure only positive integers are entered by the user
    assert input_lrg.isdigit(), "Please provide a singular positive integer"
    url = 'http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_%s.xml' % input_lrg
    r = requests.get(url, allow_redirects=True)
    open('LRG_%s.xml' % input_lrg, 'wb').write(r.content)
    fileName = open('LRG_%s.xml' % input_lrg, 'r')
    return fileName


# Alternative function to pull LRG XML from LRG website and create local XML if currently "Under Curation"
def pending_lrg_input(input_lrg):
    # Assert to ensure only positive integers are entered by the user
    assert input_lrg.isdigit(), "Please provide a singular positive integer"
    url = 'http://ftp.ebi.ac.uk/pub/databases/lrgex/pending/LRG_%s.xml' % input_lrg
    r = requests.get(url, allow_redirects=True)
    open('LRG_%s.xml' % input_lrg, 'wb').write(r.content)
    fileName = open('LRG_%s.xml' % input_lrg, 'r')
    return fileName


# Generates parsable tree of input XML file.
# Pulls whether the LRG is either currently "Published" or "Under Creation"
def tree_generation(fileName):
    try:
        tree = ET.parse(fileName)
        root = tree.getroot()
        curation = 'Curation Status: LRG Published\n'
        return tree, root, curation
    except SyntaxError:
        fileName = pending_lrg_input(input_lrg)
        tree = ET.parse(fileName)
        root = tree.getroot()
        curation = 'Curation Status: Gene Under Curation\n'
        return tree, root, curation

# Function to return exon numbers of lrg file in a list
def exon_num(root):
    exon_num_list = []
    for child in root[0]:
        if child.tag == "transcript" and child.attrib['name'] == "t1":  # Only takes the first transcirpt
            for i in child:
                if i.tag == "exon":
                    if (str(i.attrib))[-4] == "'":
                        exon_num_list.append((str(i.attrib))[-3:-2])    # Returns exons 1 - 9
                    elif (str(i.attrib))[-5] == "'":
                        exon_num_list.append((str(i.attrib))[-4:-2])    # Returns exons 10 - 99
                    elif (str(i.attrib))[-6] == "'":
                        exon_num_list.append((str(i.attrib))[-5:-2])    # Returns exons 100 - 999
    return (exon_num_list)

# Acquires the name of the gene for use in .bed file naming
def gene_name(tree):
    for gene_name in tree.findall('.//lrg_locus'):
        gene = 'LRG' + '_' + str(input_lrg) + "_" + gene_name.text
        exon_num_var = exon_num(root)
    return gene, exon_num_var


# Function to return two lists for start and end coordinates of exons respectively
def exon_coord(root):
    exon_start_list = []
    exon_end_list = []
    for child in root[0]:
        if child.tag == "transcript" and child.attrib['name'] == "t1":  # Only takes the first transcirpt
            for exon in child:
                if exon.tag == "exon":
                    exon_start_list.append(exon[0].attrib["start"])
            for exon in child:
                if exon.tag == "exon":
                    exon_end_list.append(exon[0].attrib["end"])
    return exon_start_list, exon_end_list


# Converts start_list_str and end_list_str to integer values
def list_conversion_str2int(list_a, list_b):
    output_list_a = []
    output_list_b = []
    for i in list_a[0]:
        output_list_a.append(int(i))
    for i in list_b[0]:
        output_list_b.append(int(i))
    return output_list_a, output_list_b


# Function to calculate exon lengths
def exon_len_func(lrg_start, lrg_end):
    exon_len = []
    exon_len_count = 0
    for i in lrg_start:
        exon_len.append(lrg_end[exon_len_count] - lrg_start[exon_len_count])
        exon_len_count += 1
    return exon_len


# Pulls various values from the inputted XML file
def tree_values(tree):
    for i in tree.findall('.//mapping'):
        if i.attrib["coord_system"] == "GRCh37.p13":
            chromosome = i.attrib["other_name"]  # Pulls chromosome number from XML
            gene_chr_start = int(i[0].attrib["other_start"])  # For loop to obtain start coords of gene on GRCh37.p13
            gene_chr_end = int(i[0].attrib["other_end"])  # For loop to obtain end coords of gene on GRCh37.p13
            strand = int(i[0].attrib["strand"])  # Fdentify between forward strand(1) and reverse strand(-1)
    return chromosome, gene_chr_start, gene_chr_end, strand


# Mapping LRG coords to chromosomal coordinates (FORWARD STRAND)
def strand_pos_neg(strand, lrg_start_list, lrg_end_list, gene_chr_start, gene_chr_end):
    if strand == 1:
        chr_exon_start = []
        for coord in lrg_start_list:
            chr_exon_start.append(coord + gene_chr_start - 1)
        chr_exon_end = []
        for coord in lrg_end_list:
            chr_exon_end.append(coord + gene_chr_start - 1)

    # Mapping of LRG coordinates to chromosomal locations (REVERSE STRAND)
    else:
        chr_exon_start = []
        for coord in lrg_start_list:
            chr_exon_start.append(gene_chr_end - coord + 1)
        chr_exon_end = []
        for coord in lrg_end_list:
            chr_exon_end.append(gene_chr_end - coord + 1)
    return chr_exon_start, chr_exon_end


# Pulls chromosome number from input LRG_xml
def chrom_num(chr_exon_start, chromosome):
    chr_list = []
    count = 0
    while count < len(chr_exon_start):
        chr_list.append("chr" + chromosome)
        count += 1
    return chr_list


# Creation of output BED file named by LRG # followed by gene name
def output_bed(strand, chr_list, chr_exon_start, chr_exon_end, exon_num_var, exon_len):
    date = time.strftime("File created: %d/%m/%Y  %H:%M:%S\n")  # Creates a date/time stamp for creaton of BED file
    header = "\tStart\t\tEnd\t\tExon\tLength\n"  # Headers for output text file
    if strand == 1:
        strand = "Strand: Forward (+)\n\n"
    else:
        strand = "Strand: Reverse (-)\n\n"
    # Writes generated values to the output .bed file
    with open('%s.bed' % gene, 'w+') as file_temp:
        file_temp.write(date)
        file_temp.write(curation)
        file_temp.write(strand)
        file_temp.write(header)
        for (chr_list, chr_exon_start, chr_exon_end, exon_num_var, exon_len) in zip(chr_list, chr_exon_start,
                                                                                    chr_exon_end, exon_num_var,
                                                                                    exon_len):
            file_temp.write(
                "{0}\t{1}\t{2}\t{3}\t{4}\n".format(chr_list, chr_exon_start, chr_exon_end, exon_num_var, exon_len))


# If the program is run as the __main__ program an output .bed is generated
# Allows for use of the above functions in testing processes or external programs
if __name__ == "__main__":
    input_lrg = input("Please enter LRG number: ")
    fileName = lrg_input(input_lrg)
    tree, root, curation = tree_generation(fileName)
    gene, exon_num_var = gene_name(tree)
    start_list_str, end_list_str = map(list, zip(exon_coord(root)))
    lrg_start_list, lrg_end_list = list_conversion_str2int(start_list_str, end_list_str)
    exon_len = exon_len_func(lrg_start_list, lrg_end_list)
    chromosome, gene_chr_start, gene_chr_end, strand = (tree_values(tree))
    chr_exon_start, chr_exon_end = strand_pos_neg(strand, lrg_start_list, lrg_end_list, gene_chr_start, gene_chr_end)
    chr_list = chrom_num(chr_exon_start, chromosome)
    output_bed(strand, chr_list, chr_exon_start, chr_exon_end, exon_num_var, exon_len)
