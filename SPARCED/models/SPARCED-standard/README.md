# SPARCED standard

This model is a replicate of the legacy SPARCED model.

## Antimony

The ```ant_SPARCED.txt``` file replicates the legacy Antimony SPARCED model
located in the ```Demo``` folder.
Except regarding the name of the model (addition of the prefix "ant\_"),
the two files appear similar when running:
```diff -w SPARCED.txt ant_SPARCED.txt```.
Note that the ```-w``` flag allows to ignore all white space.

## SBML

The ```sbml_SPARCED.xml``` file replicates the legacy SBML SPARCED model
located in the ```Demo``` folder.
Except that some lines were automatically written at different places in the
file, the two files appear similar when running:
```diff -w SPARCED.xml sbml_SPARCED.xml```.
Note that the ```-w``` flag allows to ignore all white space.

