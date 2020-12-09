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
    else:
      rawdata = open(filename+'.csv','r')

    data = {}
    counter = 0

    for line in rawdata:
        if line.strip():
            if counter == 0:
                counter += 1
                column_names = line.strip("\n' ''")

                if sif == True:
                  data_name = ['start_nodes','edge_features','end_nodes','receptor_switch']
                else:
                  data_name = line.split(",")

                for name in data_name:
                    data[name] = []

            else:
                line = line.strip("\n' '")
    
                if sif == True:
                  line = line.split("	")
                else:
                  line = line.split(",")
                
                c = 0

                for name in data_name:
                    element = line[c]
                    c += 1 
                    data[name].append(element)
    
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
def analysis_prep(data_network,node_names,edge_weight_dict,destination):
  counter = 0
  Receptors = []
  for switch in data_network['receptor_switch']:
    if switch == 'true':
      Receptors.append(data_network['start_nodes'][counter])
    counter += 1 
  Receptors = np.unique(Receptors)

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
      print('The shortest path from ',R,' (Upstream Receptor) to ',destination,' in a weighted Pathway Network: \n', Path,'\n with score of ',Score,'\n')
    except:
      pass

  return Receptors, pathDict, pathScore
```

## Calculating overall score 
- This will take expression from mRNA and protein expression and converts into overall score along with pathway length based score (addition)

```python
def exp_checker(Receptor_name,data_mrna,name_from_prot,prot_data,pathScore,pathDict):
  List_path_nodes = pathDict[Receptor_name]
  print('Checking the expression of each node from mRNA seq data for ', Receptor_name,'-mediated pathway')
  M2Score = 0
  PolarizeScore = 0
  prot_exp_score = 0

  for node in List_path_nodes:
    for g_name in data_mrna['mgi_symbol']:
      if node.lower() == g_name.lower():
        n = data_mrna['mgi_symbol'].index(g_name)
        expM2 = float(data_mrna['M2'][n])
        compFC = float(data_mrna['logFC(M1M2/M2)\n'][n])
        print('- Name of Node: ', node)

        if float(expM2) < 0.001:
          continue
        else:
          M2Score = M2Score + expM2
          PolarizeScore = PolarizeScore + compFC
          print('-- Expression in M2: ', expM2)
          print('-- Fold Change from M2 to M1M2: ', compFC)

    for p_name in name_from_prot:
      if node.lower() == p_name.lower():
        prot_exp = 1/prot_data[p_name]
        if prot_exp_val < 1:
          prot_exp_val = prot_exp*-1
    
        prot_exp_score += prot_exp_val
        print(p_name)
        print('-- Fold Change from Protein Expression: ',prot_exp_val)
      else:
        prot_exp_val = 0
      
  
  pathscore = pathScore[Receptor_name]
  if pathscore > 1000:
    pathscore = -1000
  else:
    pathscore = np.log10(1/pathscore)

  total = 0.1*M2Score + 10*PolarizeScore + pathscore + prot_exp_score*10
  print('- Pathway Score: ',total, 'pathway length: ',pathScore[Receptor_name])
  print('----END----')

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
def score_per_receptor(destination,data_mrna,name_from_prot,prot_data,data_network,node_names,edge_weight_dict):
  Receptors, pathDict, pathScore = analysis_prep(data_network,node_names,edge_weight_dict,destination)
  collected_score = []
  tested_receptors = []
  print(Receptors)
  for receptor_name in Receptors:
    try:
      score = exp_checker(receptor_name,data_mrna,name_from_prot,prot_data,pathScore,pathDict)
      collected_score.append(score)
      tested_receptors.append(receptor_name)
    except:
      pass

  receptor_and_score = {'Receptors': tested_receptors, 'Score': collected_score}
  print(receptor_and_score)
  return receptor_and_score
```
