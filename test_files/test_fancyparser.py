"""
PyTest script for testing the fancyparser.py LRG parser programme. To run the test the test1.xml file must be present
within the same directory. For full usage please refer to the README file.
"""

# Importing required modules
import pytest
import sys
import os
sys.path.append(os.path.abspath('../'))
import xml.etree.ElementTree as ET
import fancyparser as fp

fileName = 'test1.xml'
xml202 = 'LRG_202.xml'
input_lrg = '1'
fp.input_lrg = '1'
tree = ET.parse(fileName)
fp.root = tree.getroot()
root = tree.getroot()


# Correct naming convention, r/w mode, and file type of outputted file
def test_lrg_input():
    assert str(fp.lrg_input(input_lrg)) == "<_io.TextIOWrapper name='LRG_1.xml' mode='r' encoding='UTF-8'>"


# Correct naming convention, r/w mode, and filetype of outputted file
def test_pending_lrg_input():
    input_lrg = '14'
    assert str(fp.pending_lrg_input(input_lrg)) == "<_io.TextIOWrapper name='LRG_14.xml' mode='r' encoding='UTF-8'>"


# Assert correct ElementTree functions (tree, root) and variable (curation) are generated
def test_tree_generation():
    tree, root, curation = fp.tree_generation(fileName)
    tree = str(tree)
    root = str(root)
    curation = str(curation)
    assert tree.startswith("<xml.etree.ElementTree.ElementTree object at ")
    assert root.startswith("<Element 'lrg' at ")
    assert curation in ['Curation Status: LRG Published\n' or 'Curation Status: Gene Under Curation\n']


# Assert correct number of exons is calculated
def test_exon_num():
    tree202 = ET.parse(xml202)
    root202 = tree202.getroot()
    result = [str(i+1) for i in range(183)]
    assert fp.exon_num(root202) == result

# Assert create gene name is pulled from XML
def test_gene_name():
    gene, exon_num_var = fp.gene_name(tree)
    assert gene == 'LRG_1_COL1A1'


def test_exon_coord():
    # assert len(fp.exon_coord(root)) == 2
    assert len(fp.exon_coord(root)[0]) == len(fp.exon_coord(root)[1])
    lst = []
    for a in fp.exon_coord(root)[1] + fp.exon_coord(root)[1]:
        lst.append(int(a))
    for a in lst:
        assert str(type(a)) == "<class 'int'>"


# Conversion of affixed lists of strings to seperate lists of integers
def test_list_conversion_str2int():
    list_a = [[str(i) for i in range(20) if i %2 ==1]]
    list_b = [[str(i) for i in range(20) if i %2 ==0]]
    output_list_a, output_list_b = fp.list_conversion_str2int(list_a,list_b)
    for a in output_list_a:
        assert str(type(a))[-5:-2] == 'int'
    for b in output_list_b:
        assert str(type(b))[-5:-2] == 'int'


# Ensure subtraction of each integer in list_b from the corresponding in list_a
def test_exon_len_func():
    list_a = [i for i in range(40) if i %2 ==1]
    list_b = [i for i in range(40) if i %2 ==0]
    assert fp.exon_len_func(list_b, list_a) == [1 for i in range(len(list_a))]
    assert fp.exon_len_func(list_b, list_a) == [1 for i in range(len(list_b))]


# Testing of correct values to be returned upon parsing of the LRG XML file
def test_tree_values():
    chromosome, gene_chr_start, gene_chr_end, strand = fp.tree_values(tree)
    assert chromosome == '17'
    assert gene_chr_start == 48259457
    assert gene_chr_end == 48284000
    assert strand == -1

def test_strand_pos_neg():
    lrg_start_list = [i*3 for i in range(21)]
    lrg_end_list = [i* 3 for i in range(21)]
    gene_chr_start = 1000000
    gene_chr_end = 1000000
    # Testing positive strand mapping
    strand_pos = 1
    pos_chr_exon_start, pos_chr_exon_end = fp.strand_pos_neg(strand_pos, lrg_start_list, lrg_end_list, gene_chr_start, gene_chr_end)
    pos_count1 = 0
    pos_count2 = 0
    for i in pos_chr_exon_start:
        assert i == lrg_start_list[pos_count1] + gene_chr_start -1
        pos_count1 += 1
    for i in pos_chr_exon_end:
        assert i == lrg_end_list[pos_count2] + gene_chr_start -1
        pos_count2 += 1
    assert len(pos_chr_exon_start) == len(lrg_start_list)
    assert len(pos_chr_exon_end) == len(lrg_end_list)
    # Testing negative strand mapping
    strand_neg = -1
    neg_chr_exon_start, neg_chr_exon_end = fp.strand_pos_neg(strand_neg, lrg_start_list, lrg_end_list, gene_chr_start, gene_chr_end)
    neg_count1 = 0
    neg_count2 = 0
    for i in neg_chr_exon_start:
        assert i == gene_chr_end - lrg_start_list[neg_count1] +1
        neg_count1 += 1
    for i in neg_chr_exon_end:
        assert i == gene_chr_end - lrg_start_list[neg_count2] +1
        neg_count2 += 1
    assert len(neg_chr_exon_start) == len(lrg_start_list)
    assert len(neg_chr_exon_end) == len(lrg_end_list)



# Ensure copy of chromosome number is created for each exon
def test_chrom_num():
    chromosome = '5'
    chr_exon_start = [i+1 for i in range(30)]
    chr_list = fp.chrom_num (chr_exon_start, chromosome)
    assert len(chr_list) == len(chr_exon_start)
    for i in chr_list:
        assert i == "chr5"

"""
Remaining test functions to create:
    lrg_input                       DONE
    pending_lrg_input               DONE
    tree_generation                 DONE
    gene_name                       DONE
    exon_num                        DONE
    exon_coord                      DONE
    test_list_conversion_str2int    DONE
    exon_len_func                   DONE
    tree_values                     DONE
    strand_pos_neg                  DONE
    chrom_num                       DONE
    output_bed
"""
