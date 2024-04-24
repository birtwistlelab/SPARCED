# Project Foreman Documentation:

### Author(s): Jonah R. Huggins

Interacting with Project Foreman code is intended through the command line, where a user calls the python scripts new_test.py, unit_tests.py and createModel_unitTests.py to conduct an individualized unit test and specify user defined arguments.

## Structure and organization

Within the SPARCED model structure, specific unit test work takes place in SPARCED/unit_tests. Before conducting a unit test, the project structure is as follows:

.
├── bin
│   └── /
├── configs
│   └── /
├── docker
│   └── /
├── licenses
│   └── /
├── LICENSE.txt
├── main.nf
├── nextflow.config
├── README.md
├── requirements.txt
├── SPARCED_Brep
│   └── /
├── input_files
│   └── /
├── toolbox
│   └── /
├── TROUBLESHOOTING.md
└── unit_tests
    ├── createModel_unitTest.py
    ├── new_test.py
    └── unit_testing.py


Note: "│   └── /" signifies a directory containing other files however, for the sake of brevity, these have been removed from this tree.

 ### The new test script

The organization and structure of this work is made through the new_test.py script. Executing new_test.py creates the following organizational structure:

-   'unit_tests' directory within the SPARCED root directory. To act as a parent directory 

-   Creates a new subdirectory within unit_tests to store the appropriate files, based on the CLI name argument provided:

-   Pulls in the current input_files directory in the project's root directory

-   Creates a new directory for PEtab files, 'PEtab_files', and creates a basic template file each PEtab required file

-   Creates a new directory called scripts

-   Copies in the unit_tests.py and createModel_unitTests.py files from unit_tests  into  scripts

To execute new_test.py, simply type the following into the command line:

python3 new_test.py --name arbitrary_name

--name [STRING NOT NULL]: the name for the unit test, arbitrary_name being a user defined, arbitrary name to distinguish the unit test.

Example: unit_tests' tree structure after running python3 new_test.py --name test:

SPARCED
   │   └── /
   └── unit_tests
         ├── createModel_unitTest.py
         ├── new_test.py
         ├── test
         │   ├── input_files
         │   │   ├── Compartments.txt
         │   │   ├── excelToTSV.py
         │   │   ├── GeneReg.txt
         │   │   ├── Initializer.txt
         │   │   ├── Observables.txt
         │   │   ├── OmicsData.txt
         │   │   ├── Ratelaws.txt
         │   │   ├── ratios.txt
         │   │   ├── README.md
         │   │   ├── Species.txt
         │   │   └── StoicMat.txt
         │   ├── petab_files
         │   │   ├── conditions.tsv
         │   │   ├── measurements.tsv
         │   │   ├── observables.tsv
         │   │   ├── parameters.tsv
         │   │   └── test.yml
         │   └── scripts
         │       ├── createModel_unitTest.py
         │       └── unit_testing.py
         └── unit_testing.py


## Semantics for Condition and Observable ID's in SPARCED

-   Condition and Observable Identifiers are arbitrarily assigned, according to the PEtab format. Therefore, I use a standard naming convention for clarity's sake.

-   conditionId is assigned as follows:

-   'Species1 Model Identifier' + '_' + 'Species1 Conc.' + '_' +'Species$n Model identifier' + '_' + 'Species$n Conc.' +'_' + time


