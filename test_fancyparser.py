import pytest
import xml.etree.ElementTree as ET
import numpy as np
import time
import requests

import fancyparser as fp

# conversion of affixed lists of strings to seperate lists of integers
def test_list_conversion_str2int():
    list_a = [[str(i) for i in range(20) if i %2 ==1]]
    list_b = [[str(i) for i in range(20) if i %2 ==0]]
    output_list_a, output_list_b = fp.list_conversion_str2int(list_a,list_b)
    for a in output_list_a:
        assert str(type(a))[-5:-2] == 'int'
    for b in output_list_b:
        assert str(type(b))[-5:-2] == 'int'
    # assert len(output_list_a) + len(output_list_b) == len(list_a[0]) + len(list_b[0])
