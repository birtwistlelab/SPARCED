import amici

def import_model_from_sbml(filename, model_output_dir):
    """
    I realized the original code was pretty similar to the following one:
    https://github.com/biosimulators/Biosimulators_AMICI/blob/bf0dcb25201332e1f4a32d24873149c3800a2a4f/biosimulators_amici/core.py#L262
    Hence, as some modifications and inconsitencies were going on, I am basing
    myself on it as a model to refactor this function.
    """

    """
    filename is a string : path to SBML file
    """

    sbml_importer = amici.SbmlImporter(filename)
    sbml_model = sbml_importer.sbml




sys.path.insert(0, os.path.abspath(sbml_model_name))                        
 65     sbml_reader = libsbml.SBMLReader()                                          
 66     sbml_doc = sbml_reader.readSBML(sbml_file_name)                             
 67     sbml_model = sbml_doc.getModel()                                            
 68     sbml_importer = amici.SbmlImporter(sbml_file_name)                          
 69     constant_parameters = [params.getId() \                                     
 70                            for params in sbml_model.getListOfParameters()]      
 71     sbml_importer.sbml2amici(sbml_model_name, model_output_dir,                 
 72                              verbose=args.verbose,                              
 73                              constantParameters=constant_parameters) 
