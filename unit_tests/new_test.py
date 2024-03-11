import os
import shutil
import argparse
from typing import Optional

if __name__ == "__main__":
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description="Create directories and files.")
    parser.add_argument("--name", required=True, help="Name of the new directory")
    parser.add_argument('-C', '--custom-input-files', required=False, type=int, help='(bool) define whether custom input files are necessary', default=0)
    args = parser.parse_args()


def create_directories(name, custom_input_files: Optional[int] = None):

    #orient script within the directory
    cd = os.getcwd()
    wd = os.path.dirname(cd)

    # Create a new directory based on the provided name
    if os.path.exists(cd + '/' + name): #check if the output directory (named after the yaml file exists)
        shutil.rmtree(cd + '/' + name) # if it does, remove it ***THIS WILL DELETE EVERTHING***
        os.makedirs(cd + '/' + name) #create a new directory with the same name
        os.chdir(cd + '/' + name) 
    
    else: #if the output directory doesn't exist, create it
        os.makedirs(cd + '/' + name)
        os.chdir(cd + '/' + name)
        os.mkdir('scripts')

    # Copy create and run model files into the new directory
    shutil.copy(cd + '/src/runModel_unitTest.py', os.path.join(os.getcwd(), 'scripts/runModel_unitTest.py'))
    shutil.copy(cd + '/src/createModel_unitTest.py', os.path.join(os.getcwd(), 'scripts/createModel_unitTest.py'))
    
    #If user defines preference Copy 'input_files' directory and its contents into the new directory
    if custom_input_files == 1:
        shutil.copytree(wd + '/input_files', os.path.join(os.getcwd(), 'input_files'))

    # Create 'petab_files' directory within the new directory
    os.makedirs(os.path.join(os.getcwd(), 'petab_files'))

    # Create 4 files ('file1.txt' through 'file4.txt') within 'petab_files' directory
    petab_file_names = ['conditions', 'observables', 'parameters', 'measurements']
    for file in petab_file_names:
        filename = f'{file}.tsv'
        filepath = os.path.join(name, 'petab_files', filename)
        with open(f'petab_files/{filename}', 'a') as file:
            file.write(f'{filename}, delete on data entry')

    yaml_name = f'{name}'
    yaml_path = os.path.join(name, 'petab_files', yaml_name)
    with open(f'petab_files/{yaml_name}.yml', 'a') as file:
        file.write("""
    format_version: 1 
    parameter_file: parameters.tsv
    problems: 
      - condition_files: 
        - conditions.tsv
      measurement_files:
        - measurements.tsv
      observable_files:
        - observables.tsv
      sbml_files:
        - SPARCED.xml
    """)
        

# if __name__ == "__main__":
#     # Set up command-line argument parser
#     parser = argparse.ArgumentParser(description="Create directories and files.")
#     parser.add_argument("--name", required=True, help="Name of the new directory")

#     # Parse command-line arguments
#     args = parser.parse_args()

#     # Call the function to perform the specified tasks
create_directories(args.name, args.custom_input_files)

