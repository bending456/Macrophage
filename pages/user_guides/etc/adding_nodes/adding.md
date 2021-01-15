---
title: How to add new members     
nav_order: 2
parent: Etc.
grand_parent: User Guides
---

# Adding new members to the network and assocaited analyses. 
When adding a new ligand/receptor mediated pathway, you will need to update the following files at a minimum:
networkXXX.sif 
- add the ligand (GENE1) and receptor (GENE2), as well as any downfield node that the latter controls (GENE3)
GENE1 up-regulates GENE2
GENE2 up-regulates GENE3

receptorlist.txt 
- add receptor GENE2

ligandlist.yaml 
- add receptor/ligand pair
GENE2: [GENE1]

