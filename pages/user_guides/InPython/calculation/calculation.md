---
title: Calculation in Detail
nav_order: 2
parent: Colaboratory
grand_parent: User Guides
---

# What is in **utilities.py**?
- This file will be downloaded along with data file in /home/pathwayanalysis once the first cell of script is executed. 

## Loading packages required for the calculation/analysis
```python
import numpy as np
from networkx.drawing.nx_agraph import graphviz_layout
from graphviz import Digraph
import networkx as nx 
from IPython.display import Image
```

## File reader for mRNA sequence data and Cytoscape pathway network files
- mRNA sequence data: CSV file 
- pathway network data: SIF file 
```python
def reader(filename,
           sif=False):

    if sif == True:

      rawdata = open(filename+'.sif','r')
      data = {}
      network_info = ['start_nodes','edge_features','end_nodes']
      counter = 0

      for line in rawdata:
        if line.strip():
          line = line.strip("\n' '")
          line = line.split("	")

          if counter == 0:
            counter = counter + 1 
            for name in network_info:
              data[name]=[]

          for i in np.arange(len(network_info)):
            data[network_info[i]].append(line[i])

      rawdata.close()

    else:

      rawdata = open(filename+'.csv','r')
      data = {}
      counter = 0

      for line in rawdata:
        if line.strip():
          if counter == 0:
            counter = counter + 1 
            column_names = line.strip("\n' '")
            data_name = line.split(",")
            # Creating the dictionary data structure
            for name in data_name:
              data[name]=[]
          else:
            line = line.strip("\n' '")
            line = line.split(",")
            # Storing the data in the dictionary 
            for i in np.arange(4):
              element = line[i]
              data[data_name[i]].append(element)

    
      rawdata.close()
    return data
```

## File reader for List of Receptors
``` python
def textreader(filename):
    rawdata = open(filename+'.txt','r')

    data = []
    counter = 0

    for line in rawdata:
        if line.strip():
          line = line.strip("\n' '")
          data.append(str(line))
    
    rawdata.close()
    return data
```

## File reader for protein expression data
- CSV file 
```python
def protReader(filename):
  rawdata = open(filename+'.csv','r')
  data = {}
  name = []
  counter = 0 

  for line in rawdata:
    if line.strip():
      if counter == 0:
        counter = counter + 1
        continue 
      else:
        line = line.strip("\n' '")
        line = line.split(",")
        # Storing the data in the dictionary 
        if line[1] == "N/A":
          continue
        else:
          data[line[1]] = float(line[2])
          name.append(line[1])

  return data, name
```

## Assigning edge weight based on the name of edge 
- In this specific example, we use "up-regulation" and "down-regulation" only. However, if the future process can elaborate each edge, we recommend to sophisticate this process. Please, feel free to contact us! 

```python
def weighingEdges(edge_names,down_score,up_score):
    edge_weight_dict = {}
    for name in edge_names:
        if name == 'down-regulates':
            w = [down_score,'D']
            edge_weight_dict[name] = w
        elif name == 'up-regulates':
            w = [up_score,'U']
            edge_weight_dict[name] = w

    return edge_weight_dict
```

## Generating network figure via NetworkX 
- The png file will be stored in your current working directory. 
- In this particular example, the current working directory is assigned to be /content/pathwayanalysis
- The file name is 'network_figure.png'

```python
def NetworkFigure(node_names,data_network,edge_weight_dict):
  G2 = Digraph('unix', filename='fullpicture',format='png',
            node_attr={'shape':'box',
                       'color': 'blue:purple', 
                       'style': 'filled',
                       'fontcolor':'white'},
            edge_attr={'color': 'red'})
  
  for node in node_names:
    G2.node(node)

  for n in np.arange(len(data_network['start_nodes'])):
    weight = edge_weight_dict[data_network['edge_features'][n]]
    G2.edge(data_network['start_nodes'][n],data_network['end_nodes'][n],label=str(weight))

  G2.view('network_figure')

  return
```

## Analysis prep function 1 
- This function spits out:
    1. List of receptors in pathway (list of target nodes as starting nodes)
    2. List of nodes in the pathway led by each receptor 
    3. List of score calculated from path length based on edge weight

```python
def analysis_prep(data_network,
                  receptorlist,
                  node_names,
                  edge_weight_dict,
                  destination):

  counter = 0
  Receptors = receptorlist

  G = nx.DiGraph()
  for node in node_names:
    G.add_node(node)
  
  G.nodes()

  for n in np.arange(len(data_network['start_nodes'])):
    weight = edge_weight_dict[data_network['edge_features'][n]][0]
    label_n = edge_weight_dict[data_network['edge_features'][n]][1]
    G.add_edge(data_network['start_nodes'][n],data_network['end_nodes'][n],weight=weight,label=label_n)   
  G.edges()

  Outcome = destination 
  pathDict = {}
  pathScore = {}

  for R in Receptors:
    try:
      Path = nx.dijkstra_path(G,R,Outcome)
      pathDict[R] = Path
      Score = nx.dijkstra_path_length(G,R,Outcome)
      pathScore[R] = Score
      #print('The shortest path from ',R,' (Upstream Receptor) to ',destination,' in a weighted Pathway Network: \n', Path,'\n with score of ',Score,'\n')
    except:
      pass

  return Receptors, pathDict, pathScore
```

## Calculating overall score 
- This will take expression from mRNA and protein expression and converts into overall score along with pathway length based score (addition)

```python
def exp_checker(Receptor_name,
                data_mrna,
                name_from_prot,
                prot_data,
                pathScore,
                pathDict):

  List_path_nodes = pathDict[Receptor_name]
  #print('Checking the expression of each node from mRNA seq data for ', Receptor_name,'-mediated pathway')
  #print(pathDict[Receptor_name])
  M2Score = 0
  PolarizeScore = 0
  prot_exp_score = 0

  for node in List_path_nodes:
    for g_name in data_mrna['mgi_symbol']:
      if node.lower() == g_name.lower():
        n = data_mrna['mgi_symbol'].index(g_name)
        expM2 = float(data_mrna['M2'][n])
        expM1 = float(data_mrna['M1'][n])
        expM1M2 = float(data_mrna['M1M2'][n])
        compFC = math.log2(expM1M2/expM2)
        #print('- Name of Node: ', node)

        if float(expM2) < 0.001:
          continue
        else:
          M2Score = M2Score + expM2
          PolarizeScore = PolarizeScore + compFC
          #print('-- Expression in M2: ', expM2)
          #print('-- Fold Change from M2 to M1M2: ', compFC)

    for p_name in name_from_prot:
      if node.lower() == p_name.lower():
        prot_exp = 1/prot_data[p_name]
        if prot_exp_val < 1:
          prot_exp_val = prot_exp*-1
    
        prot_exp_score += prot_exp_val
        #print(p_name)
        #print('-- Fold Change from Protein Expression: ',prot_exp_val)
      else:
        prot_exp_val = 0
      
  
  pathscore = pathScore[Receptor_name]
  if pathscore > 1000:
    pathscore = -1000
  else:
    pathscore = np.log10(1/pathscore)

  total = 0.1*M2Score + 10*PolarizeScore + pathscore + prot_exp_score*10
  #print('- Pathway Score: ',total, 'pathway length: ',pathScore[Receptor_name])
  #print('----END----')

  return total
```

- The following is the way converting expression into score
```python
total = 0.1*M2Score + 10*PolarizeScore + pathscore + prot_exp_score*10
```
- *pathscore* is either quite small (< 100) or big (> 1000) number. Small number indicates more favorable path than big score pathway because 1000 score is added to the pathway score if path is inhibitory. Therefore, we convert the score to -1000 if the original score is bigger than 1000. 
- *M2Score* indicates the expression of M1 specific gene in M2 under the assumption that its abundance will be addition for M1 polarization. 
- *PolarizeScore* indicates the expression change of M1 specific gene due to M1 MEV treatment (M1M2). 
- *prot_exp_score* indicates the expression upregulated from M2 to M1. Only selective proteins listed as nodes in selected-receptor mediated pathway. 

## Generating the final outcome. 
```python
def score_per_receptor(destination,
                       data_mrna,
                       receptor_list,
                       name_from_prot,
                       prot_data,
                       data_network,
                       node_names,
                       edge_weight_dict):

  Receptors, pathDict, pathScore = analysis_prep(data_network,
                                                 receptor_list,
                                                 node_names,
                                                 edge_weight_dict,
                                                 destination)
  collected_score = []
  tested_receptors = []

  #print(Receptors)
  for receptor_name in Receptors:
    try:
      score = exp_checker(receptor_name,
                          data_mrna,
                          name_from_prot,
                          prot_data,
                          pathScore,
                          pathDict)

      collected_score.append(score)
      tested_receptors.append(receptor_name)
    except:
      pass

  receptor_and_score = {'Receptors': tested_receptors, 'Score': collected_score}
  #print(receptor_and_score)
  return receptor_and_score, pathDict
```

## Checking mRNA sequence and protein expression data against specific node
```python
def expression_search(name_node,
                      data_mrna,
                      prot_data,
                      name_from_prot):
  
  for name_prot in name_from_prot:
    if name_node.lower() == name_prot.lower():
      element_prot = prot_data[name_prot]
  
  try:
    if element_prot:
      prot_exp = element_prot
  except:
    prot_exp = 'N/A'
  
  for name_mrna in data_mrna['mgi_symbol']:
    if name_node.lower() == name_mrna.lower():
      n = data_mrna['mgi_symbol'].index(name_mrna)
      expM2 = float(data_mrna['M2'][n])
      expM1 = float(data_mrna['M1'][n])
      expM1M2 = float(data_mrna['M1M2'][n])
      compFC = math.log2(expM1M2/expM2)
  
  try:
    if n:
      M1 = expM1
      M2 = expM2 
      M1M2 = expM1M2 
      logFCofM1M2overM2 = compFC 
  except:
    M1 = 'N/A'
    M2 = 'N/A'
    M1M2 = 'N/A'
    logFCofM1M2overM2 = 'N/A'

  return prot_exp, M1, M2, M1M2, logFCofM1M2overM2
  ```


# What's in ***pathwaySearch.py***?
```python
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout
from graphviz import Digraph
import utilities as util
import yaml 

def runItAll(
              mRNA_inputfilename = 'mRNAdata', 
              Prot_inputfilename = 'Ab_Chris',
              network_inputfilename = 'network',
              receptor_list = 'receptorlist',
              List_file_name = 'test2',
              destination = 'M1_polarization',
              unit_test_receptor = None
              ): 
  #####################################################
  ### Extracting Gene names from mRNA sequence data ###
  #####################################################
  data_mrna = util.reader(mRNA_inputfilename)

  #########################################
  ### Extracting protein names and data ###
  #########################################
  data_prot, prot_name = util.protReader(Prot_inputfilename)

  ########################################################################
  ### Extracting Node names from pathway file generated from Cytoscape ###
  ########################################################################
  data_network = util.reader(network_inputfilename,sif=True) 

  #########################################################
  ### Extracting Receptor node names from Receptor List ###
  #########################################################
  receptor_list = util.textreader(receptor_list) 

  ##############################################################################
  ### Sorting nodes and edges                                                ###
  ###           and add the weight to each edge based on its unique feature  ###
  ##############################################################################
  node_names = np.unique(data_network['start_nodes']+data_network['end_nodes']) ## get unique node name
  edge_unique_feature = np.unique(data_network['edge_features']) ## get unique edge feature

  edge_weight_dict = util.weighingEdges(edge_unique_feature,1000,1)

  ######################################################################
  ###       Generating Network Image based on pathway file           ###
  ###  This step is simply for the better visualization than default ###
  ######################################################################
  '''
  File will be generated in /home/pathwayanalysis directory
  under the name "network_figure.png"
  '''
  util.NetworkFigure(node_names,data_network,edge_weight_dict)

  ###########################################
  ###       Performing the analysis       ###
  ###########################################
  receptor_and_score, pathDict = util.score_per_receptor(destination,
                                                         data_mrna,
                                                         receptor_list,
                                                         prot_name,
                                                         data_prot,
                                                         data_network,
                                                         node_names,
                                                         edge_weight_dict)

  with open(List_file_name+'.yaml') as file:
    ligand_list = yaml.load(file)

  receptor_specifics = {}
  ligand_specifics = {}

  for R in receptor_list:
    receptor_specifics[R] = {}
    ligand_specifics[R] = {}

    try:
      list_node = pathDict[R]
      for node in list_node:
        receptor_specifics[R][node] = {}
        prot_exp, M1, M2, M1M2, logFCofM1M2overM2 = util.expression_search(node,
                                                                           data_mrna,
                                                                           data_prot,
                                                                           prot_name)

        receptor_specifics[R][node]['mAb'] = prot_exp
        receptor_specifics[R][node]['M1_mRNA'] = M1
        receptor_specifics[R][node]['M2_mRNA'] = M2
        receptor_specifics[R][node]['M1M2_mRNA'] = M1M2
        receptor_specifics[R][node]['logFC(M1M2/M2)'] = logFCofM1M2overM2
       
    except:
      pass

    try: 
      list_ligand = ligand_list[R]    
      for ligand in list_ligand:
        ligand_specifics[R][ligand] = {}
        prot_exp, M1, M2, M1M2, logFCofM1M2overM2 = util.expression_search(ligand,
                                                                           data_mrna,
                                                                           data_prot,
                                                                           prot_name)

        ligand_specifics[R][ligand]['mAb'] = prot_exp
        ligand_specifics[R][ligand]['M1_mRNA'] = M1
        ligand_specifics[R][ligand]['M2_mRNA'] = M2
        ligand_specifics[R][ligand]['M1M2_mRNA'] = M1M2
        ligand_specifics[R][ligand]['logFC(M1M2/M2)'] = logFCofM1M2overM2

    except:
      pass

  if unit_test_receptor:
    Receptors, pathDict, pathScore = util.analysis_prep(data_network,
                                                        receptor_list,
                                                        node_names,
                                                        edge_weight_dict,
                                                        destination)

    score = util.exp_checker(unit_test_receptor,
                             data_mrna,
                             prot_name,
                             data_prot,
                             pathScore,
                             pathDict)
  else:
    score = 'This is not unit test'

  import pandas as pd 
  df = pd.DataFrame(receptor_and_score)

  return df, receptor_specifics, ligand_specifics, receptor_list, score

#!/usr/bin/env python
import sys
##################################
#
# Revisions
#       10.08.10 inception
#
##################################


#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -run" % (scriptName)
  msg+="""
  
 
Notes:

"""
  return msg

#
# MAIN routine executed when launching this script from command line 
#
if __name__ == "__main__":
  import sys
  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  #fileIn= sys.argv[1]
  #if(len(sys.argv)==3):
  #  1
  #  #print "arg"

  # Loops over each argument in the command line 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-run"):
         runItAll()
         quit()

  raise RuntimeError("Arguments not understood")
  ```