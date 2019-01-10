# Ensure a local copy of "LRG_1.xml" exists before runinig pytest to ensure correct functionality of this script

import xml.etree.ElementTree as ET
import numpy as np
import time
import requests

import fancyparser as fp

# Correct naming convention, r/w mode, and filetype of outputted file
def test_lrg_input():
    input_lrg = '1'
    assert str(fp.lrg_input(input_lrg)) == "<_io.TextIOWrapper name='LRG_1.xml' mode='r' encoding='UTF-8'>"

# Correct naming convention, r/w mode, and filetype of outputted file
def test_pending_lrg_input():
    input_lrg = '14'
    assert str(fp.pending_lrg_input(input_lrg)) == "<_io.TextIOWrapper name='LRG_14.xml' mode='r' encoding='UTF-8'>"

fileName = 'LRG_1.xml'
tree = ET.parse(fileName)
root = tree.getroot()

# conversion of affixed lists of strings to seperate lists of integers
def test_list_conversion_str2int():
    list_a = [[str(i) for i in range(20) if i %2 ==1]]
    list_b = [[str(i) for i in range(20) if i %2 ==0]]
    output_list_a, output_list_b = fp.list_conversion_str2int(list_a,list_b)
    for a in output_list_a:
        assert str(type(a))[-5:-2] == 'int'
    for b in output_list_b:
        assert str(type(b))[-5:-2] == 'int'

# testing of correct values to be returned upon parsing of the LRG XML file
def test_tree_values():
    chromosome, gene_chr_start, gene_chr_end, strand = fp.tree_values(tree)
    assert chromosome == '17'
    assert gene_chr_start == 48259457
    assert gene_chr_end == 48284000
    assert strand == -1

"""
Remaining test functions to create:
    tree_generation
    gene_name
    exon_num
    exon_coord
    exon_len_func
    strand_pos_neg
    chrom_num
    output_bed
"""
