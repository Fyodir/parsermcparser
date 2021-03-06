## XML Parser for LRG files

An XML parser in Python for producing a .bed file from an LRG_xml file obtained from the LRG-Sequence website
"https://www.lrg-sequence.org".

*n.b. This program is designed at this time to only work with the GRCh37.p13 Human Reference Genome build*

## Authors
- **Andrew Smith**
- **Jethro Rainford**

<br/>

Development start date: 27th November 2018.

Available at: https://github.com/Fyodir/parsermcparser .

## Prerequisites


- Python (3.6.x)
- Modules
    - elementTree (from XML library) (Standard Library)
    - time (Standard Library)
    - requests (2.21.0)
    - pytest (for testing) (4.1.0)
- Working internet connection

---

## Usage

First run  `pip install --upgrade -r requirements.txt` to ensure that all correct modules are installed on your system
prior to use of the `fancyparser.py` program. If the required modules are not present, the prior command will install
them onto your system.

### fancyparser.py

Utilises a user input to return an LRG.xml file from the LRG website for use in generating
a .bed file.

##### Input:

Use input of desired LRG number.

##### Output:

.xml file of queried gene titled "LRG_(LRG number here).xml" (i.e. `LRG_1.xml`).

.bed file of queried gene titled "LRG_(LRG number)(HGNC nomenclature).bed" (i.e. `LRG_1_COL1A1.bed`).

##### Instructions for use:

- Navigate to the ```/parsermcparser``` directory on the bash terminal.
- Run the cmd line `````python3 fancyparser.py`````.
- You will be prompted with  ```Please enter LRG number: ``` request.
- Enter LRG number of desired query gene (LRG numbers for HGNC names may be found at: 
  https://www.lrg-sequence.org/search/?query=* ).
- The .xml file for the requested LRG and generated .bed file will be returned to the current directory, and may be 
  viewed with a standard text editor.

##### Example usage:
```
python3 fancyparser.py [enter]
Please enter LRG number: 39 [enter]
```
<br/>

---

## Testing

- The correct functioning of the parser was determined by the use of assert statements and error generation throughout
  the script.
- A number of modified, faulty, and corectly formatted test XML files were used to ensure functioning of test features
  - test1.xml, LRG1.xml, LRG_14.xml, LRG_202.xml.
- Assert statements have been included to reject incorrect file types, multiple input LRG numbers, incorrectly formatted
  / damaged XML files.
- The standard PyTest testing framework was also used to test the whole script for correct functioning (see below).

### PyTest
The module "PyTest" was used to test the functionality of the code and to ensure that each function works as intended.

The test files are located within a separate ```test_files``` directory.

#### Usage:

- Ensure you are within the ```test_files``` directory.
- Ensure both ```test_fancyparser.py``` and ```test1.xml``` are present (test.xml is an unmodified copy of the 
  LRG_1.xml, renamed for clarity).
- To invoke testing execute the command ```pytest -rpf```.
- This will search within the directory for the test script and run on the ```test1.xml``` file.
- The flags ```-r``` (short summary), ```-p``` (passed) and ```-f``` (fail) are recommended to give a basic visual
  output of the passed/failed tests, other variables may be found in the relevant ```man``` page.
