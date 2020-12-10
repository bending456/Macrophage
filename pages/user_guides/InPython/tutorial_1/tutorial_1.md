---
title: Colab Version Script
nav_order: 1
parent: Colaboratory
grand_parent: User Guides
---

# What's in this script?
![This is the main page of Colab version script](main.png)
## Getting started
- [Colab version script](https://colab.research.google.com/drive/1PSOWYYL02cIj1OOGdbbzifEe-YYkp6cL?usp=sharing): 
- Make sure you copy this file into your own Google Drive
  - Click "File" in the main tap
  - Click "Save a copy in Drive"

## Loading Packages
- Execute the following script in ipynb or Colab script by either click the play button or press [Shift] + [Enter]

```python
import os
os.system('git clone https://bending456@bitbucket.org/pkhlab/pathwayanalysis.git')
os.chdir('/content/pathwayanalysis') ## <------ This is your working directory so click the folder symbol to navigate it

os.system('apt-get install graphviz libgraphviz-dev pkg-config')
!pip3 install networkx
!pip3 install pygraphviz
!pip3 install graphviz

import networkx as nx 
import numpy as np 
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout
from graphviz import Digraph
import utilities as util
import pathwaySearch as ps
import pandas as pd 

%load_ext autoreload
%autoreload 2
```

## Analysis
### User Guide 
- The first outcome as you execute the following cell is the list of receptor names in the network library and corresponding score that the script processed  

```python
dataSets, receptor_specifics, receptors  = ps.runItAll(
                                                mRNA_inputfilename = 'mRNAdata', 
                                                Prot_inputfilename = 'Ab_Chris',
                                                network_inputfilename = 'network',
                                                receptor_list = 'receptorlist',
                                                destination = 'M1_polarization')
dataSets.sort_values('Score',ascending=False)
```
- ***[Note]*** Check out the data structure at [Here](/pages/user_guides/etc/data_structure/data_structure.html)
- mRNA_inputfilename: the name of file containing mRNA sequence data
- Prot_inputfilename: the name of file containing Protein Expression data  
- network_inputfilename: the name of network file (sif)
- receptor_list: the name of file containing the list of receptors in network 
- destination: the target node from the receptor 

### Outcome 1 
![This is the first outcome](py1.png)

### Checking out the raw data 
- The following cell will print out the table based on a name of receptor you desire (Receptor Names listed in the table above). 

### Column 1: Name of data 
- mAb: protein expression data (M2/M1) 
- M1_mRNA: mRNA expression data in M1 
- M2_mRNA: mRNA expression data in M2  
- M1M2_mRNA: mRNA expression data in M1M2 
- logFC(M1M2/M2): Fold Change (base 2) from M2 to M1M2 (MEV treated M2)

#### Raw 1: Name of nodes

To obtain the new table from the following example, type a name of receptor 

```python

###############################################
# ------------------------------------------- #
#                                             #
# Receptor_Name = 'Name_of_Receptor_you_want' #
#                                             #
# --------------------------------------------#
###############################################
Receptor_Name = 'IL6R'
sorted = pd.DataFrame(receptor_specifics[Receptor_Name])
sorted
```

### Outcome 2
![This is the second outcome](py2.png)

### Locating Visulaized Network image 
- In the directory of ***/content/pathwayanalysis*** (check out the first image posted at top), you should be able to locate ***network_figure.png***

![This is network](network_figure.png)