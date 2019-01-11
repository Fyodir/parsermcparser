# Ensure a local copy of "LRG_1.xml" exists before runinig pytest to ensure correct functionality of this script

import xml.etree.ElementTree as ET
import fancyparser as fp

fileName = 'LRG_1.xml'
input_lrg = '1'

# Correct naming convention, r/w mode, and filetype of outputted file
def test_lrg_input():
    assert str(fp.lrg_input(input_lrg)) == "<_io.TextIOWrapper name='LRG_1.xml' mode='r' encoding='UTF-8'>"

# Correct naming convention, r/w mode, and filetype of outputted file
def test_pending_lrg_input():
    input_lrg = '14'
    assert str(fp.pending_lrg_input(input_lrg)) == "<_io.TextIOWrapper name='LRG_14.xml' mode='r' encoding='UTF-8'>"

def test_tree_generation():
    tree, root, curation = fp.tree_generation(fileName)
    tree = str(tree)
    root = str(root)
    curation = str(curation)
    assert tree.startswith("<xml.etree.ElementTree.ElementTree object at ")
    assert root.startswith("<Element 'lrg' at ")
    assert curation == "Curation Status: LRG Published" or "Curation Status: Gene Under Curation"

tree = ET.parse(fileName)
root = tree.getroot()

# def test_gene_name():
#     gene, exon_num_var = fp.gene_name(tree)
#     assert gene == "LRG_1_COL1A1"

def test_exon_num():
    result = [str(i+1) for i in range(51)]
    assert fp.exon_num(root) == result


def test_exon_coord():
    # assert len(fp.exon_coord(root)) == 2
    assert len(fp.exon_coord(root)[0]) == len(fp.exon_coord(root)[1])
    lst = []
    for a in fp.exon_coord(root)[1] + fp.exon_coord(root)[1]:
        lst.append(int(a))
    for a in lst:
        assert str(type(a)) == "<class 'int'>"



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
    lrg_input                       DONE
    pending_lrg_input               DONE
    tree_generation                 DONE
    gene_name
    exon_num                        DONE
    exon_coord                      DONE
    test_list_conversion_str2int    DONE
    exon_len_func
    tree_values                     DONE
    strand_pos_neg
    chrom_num
    output_bed
"""
