# Algorithm for traversing and scoring polarization pathways in a cytoscape-curated network. 

See [link](https://bending456.github.io/Macrophage/) for instructions. 



# Notes
When manually-editing a .sif file, make sure that tabs are used (e.g.   GENE1<tab>up-regulates<tab>GENE2) 

When adding a new ligand/receptor mediated pathway, you will need to update the following files at a minimum:
networkXXX.sif - add the ligand (GENE1) and receptor (GENE2), as well as any downfield node that the latter controls (GENE3)
GENE1 up-regulates GENE2
GENE2 up-regulates GENE3
receptorlist.txt - add receptor GENE2
ligandlist.yaml - add receptor/ligand pair
GENE2: [GENE1]
