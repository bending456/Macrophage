---
title: How to add new members     
nav_order: 2
parent: Etc.
grand_parent: User Guides
---

# What are the given data?
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


# Notes
- Use Alpha numerical characters only! (PLCB, not PLCβ)
- Only use GENE name, not protein name (G proteins are named as complexes right now, e.g. Gai, but need to figure this out later). 
- When adding nodes, one of the intermediaries needs to connect to a preexisting part of the network 
- search for names that are already in the list (such as CXCR4) to ensure you're not adding duplicate entries
- Be very careful about the order (CXCR4 upreg Gai, not Gai upreg CXCR4)
- Double check pathways you've added in cytoscape to make sure they're consistent with KEGG, etc
- Make sure all entries are connected to a preexisting node in the network


