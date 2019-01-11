## XML Parser for LRG files

An XML parser in Python for producing a .bed file from an LRG_xml file obtained from the LRG-Sequence website "https://www.lrg-sequence.org/index.html"

*n.b. This program is designed at this time to only work with the GRCh37.p13 Human Reference Genome build*

## Authors
- **Andrew Smith**
- **Jethro Rainford**

<br/>

Development start date: 27th November 2018.

## Prerequisites


- Python (3.x.x)
- Modules
    - elementTree (from XML library) (Standard Library)
    - time (Standard Library)
    - numpy (1.15.4)
    - requests (2.21.0)
    - pytest (for testing) (4.1.0)
- Working internet connection

---

## Usage

First run the `pip install -r requirements.txt` to ensure that all correct modules are installed on your system prior to use of the `fancyparser.py` program. If the required modules are not present, the prior command will install them onto your system

### fancyparser.py

Utilises a user input to extract an LRG.xml file from the LRG website for use in generating
a BED file, rather than requiring the xml file to already be located within the "parsermcparser" directory

##### Input:

Use input of desired LRG number

##### Output:

.xml file of queried gene titled "LRG_(LRG number here).xml"
.bed file of queried gene titled "LRG_(LRG number here)_(HGNC nomenclature).bed"

##### Instructions for use:

- Navigate to the ```/parsermcparser``` directory on the bash terminal
- Run the cmd line `````python3 fancyparser.py`````
- You will be prompted with the ```Please enter LRG number: ``` request
- Enter LRG number of desired query gene
- View the outputted BED file via either a text editor (ie gedit, nano) or call the name of the file using bash cmd "cat"

##### Example Input Sequence:
```
python3 fancyparser.py [enter]
Please enter LRG number: 39 [enter]
```
<br/>

---

### parserv1.py

*(Retired original working version of the LRG XML parsing program)*

Target XML file must be located within the ```/parsermcparser``` file directory

##### Input:

pass input xml file to program as string input

##### Output:

.bed file of queried gene titled ```LRG_(LRG number here)_(HGNC nomenclature).bed```

##### Instructions:
- Ensure target LRG_xml file is located within the ```/parsermcparser``` directory
- Naviagte to the ```/parsermcparser``` directory on the bash terminal
- Run the cmd line ```python3 parserv1 (insert target xml file) [enter]```
- A bed file will be created in the ```/parsermcparser``` directory titled: ```LRG_(LRG number here)_(HGNC nomenclature).bed```

##### Example Input Sequence:
```
python3 parserv1 ./LRG_391.xml [enter]
```
<br/>

---

## Testing

- The correct functioning of the parser was determined by the use of assert statements and error generation throughout the code.
- Modified faulty test XML files were used to ensure functioning of test features
- Assert statements have been included to reject incorrect file types, multiple input LRG numbers, incorrectly formatted / damaged XML files

### pytest
The module "pytest" was used to test the functionality of the code and to ensure that each function works as intended.

A second `.py` file ```test_fancyparser.py``` exists within this package that when used in conjunction with the pytest module ensures the correct functionality of the ```fancyparser.py``` script

#### Instructions:
- Ensure a copy of the ```LRG_1.xml``` file exists locally, directly  within the ```/parsermcparser``` directory
- Run the cmd ```pytest``` from the bash cmd line when positioned within the ```/parsermcparser``` directory
