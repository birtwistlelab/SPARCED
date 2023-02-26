# SPARCED Input Data

## Formatting

The data in these files are all tab-delimited between fields. 

## Excel

If you have an excel file formatted to have the same number of columns as its associated `.txt` file in this folder, it can be converted using `excelToTSV.py` provided here.



## Input files: 

These are a set of tab-delimited text files describing the various components needed to build the model. In order to perform modifications or expansions of the model, the user may edit these input files with adequate specifications of relevant information while keeping the structure intact. 

	OmicsData: This file specifies genomic, transcriptomic and proteomic information about the cell context as well as several important parameters for the gene expression module. Each row corresponds to one gene and the columns contain different data types. The information in the columns in the order of their appearance include gene names (HGNC identifiers), gene copy number, mRNA copy number per cell (mpc), gene activation rate, gene inactivation rate, constitutive transcription rate constants, maximal transcription rate constants, mRNA degradation rate, protein copy number (mpc), protein half life and translation rate constants.  

	Species: This file contains information about the species in the deterministic protein interactions module. The rows represent species including proteins, protein complexes, post-translationally modified species and mRNAs. The information in the columns includes species names, species home compartment, initial species concentration in nM and ENSEMBL gene identifiers mapped to each species. 

	Ratelaws: This file contains information about the reactions, their rate laws and corresponding rate constant values. Each row refers to a reaction in the model. The columns contain the reaction name, a description of the rate law (optional, intended for reactions modelled without the mass action assumption), and value of each rate constant appearing in the rate law.   

	StoicMat: It is the stoichiometric matrix for the deterministic reaction network. It is a s × r matrix whereby s and r are the number of species and reactions. Any element at the index [i,j] of the matrix represents the stoichiometric coefficient of the species i  with respect to reaction j.    

	GeneReg: This file describes the parameters for the species that are subject to transcriptional regulation. Here the rows correspond to individual genes and the columns are associated with several transcriptional activators or repressors that may regulate the protein species mapped to the said genes. The gene names are expressed in HGNC format in the same manner as the OmicsData input file. A cell may have zero or non-zero values. The zero values indicate no effect of the corresponding regulator to the gene while the non-zero values have the format “A;B”. Here, A is the hill coefficient and B is the half-maximal concentration of the regulated species. Positive and negative value of A represent activator and repressor respectively. 

	Compartments: This file contains information about the compartments of the model, namely, the compartment name, its volume in liters, and the corresponding GO-term 

	Observables: Some of the important species of the SPARCED model are proteins in various forms, i.e., nascent proteins, various post-translational modifications, and protein complexes. In order to be able to compare the simulation outputs to proteomics data, unique proteins across various species have to be summed up. There are 102 functionally unique proteins in total in SPARCED model, which we refer to as “protein conglomerates”. The “Observables” input file describes how model species are to be summed up to calculate protein conglomerates. In this case, each “observable” corresponds to the compartmental-volume-corrected summation of across all formats of protein species for a protein conglomerate. The rows correspond to model species and the columns correspond to observables, of which there are 102. The entries represent the number of molecules of each species to be summed in each observable.  

	Initializer (optional): This file contains information used for model initialization. Species concentrations (columns 1–2), mRNA level adjustments (columns 3–4), parameter values (columns 5–7), observables to exclude from translation rate adjustments (column 8), and single parameter scan range (columns 9–11) are populated for each step of initialization. 
