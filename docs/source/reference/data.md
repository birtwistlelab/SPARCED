# Data

Input data files are written in TSV (tab separated value) format.
A conversion script from Excel to TSV is available in the
[general utilities](utils.rst).

## Core

The core subfolder is dedicated to input files necessary for model creation
(core SPARCED).

### Compartments

The columns are corresponding to:

  - Compartment's name
  - Compartment's volume
  - Annotation in GO terms

### Observables

### Ratelaws

The columns are corresponding to:

  - Reaction's name
  - Compartment's correction
  - A single coefficient for mass action laws, or the full reaction for more
  elaborated laws
  - Values of the coefficients you passed in the law, by order of appearance in
  the equation

### Species

The columns are corresponding to:
  - Specie's name
  - Species location, which must be a valid compartment
  - Initial concentration of the specie
  - Annotation in ENSEMBL
  - Annotation in HGNC

### StoicMat

**WARNING:** This file is _huge_, do not attempt to open with anything else but
Excel, or be prepared for a crash. You have been warned.

## Simulation

The simulation subfolder is dedicated to input files necessary to run a SPARCED
model simulation.

### ligands

The columns are corresponding to:

  - Ligand's human-readable name
  - Ligand's name, which must be a valid specie
  - Ligand's initial concentration

### GeneReg

Don't touch this.

### OmicsData

Don't touch this.

