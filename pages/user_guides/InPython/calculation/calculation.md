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
          line = line.split()

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

## Assigning edge weight based on the expression level reported in mRNA sequence data set

```python
def prepare_process(networkFilename, mRNAseqFilename, protFilename = None):
    networkdata = reader(networkFilename, sif=True)
    rnadata = reader(mRNAseqFilename)

    if protFilename:
        protdata, protname = protReader(protFilename) # this will be fixed but for now, nothing will happen 

    num_of_edges = len(networkdata['start_nodes'])
    unique_start = np.unique(networkdata['start_nodes'])
    unique_end = np.unique(networkdata['end_nodes'])
    unique_all = np.unique(networkdata['start_nodes'] + networkdata['end_nodes'])

    M2_mRNA_Node = {}
    M2_mRNA_Node_arr = []  
    M1M2_mRNA_Node = {}
    M1M2_mRNA_Node_arr = []
    M1_mRNA_Node = {}
    M1_mRNA_Node_arr = []

    # FC from M2 to M1M2 (specific case for Jill's project )
    logFC_mRNA_Node_M2toM12 = {}
    logFC_mRNA_Node_arr_M2toM12 = []
    # FC from M12 to M2 (specific case for Jill's project )
    logFC_mRNA_Node_M12toM2 = {}
    logFC_mRNA_Node_arr_M12toM2 = []
    # FC from M2 to M1 (M1 polarization from M2, general case)
    logFC_mRNA_Node_M2toM1 = {}
    logFC_mRNA_Node_arr_M2toM1 = []
    # FC from M1 to M2 (M2 polarization from M1, general case)
    logFC_mRNA_Node_M1toM2 = {}
    logFC_mRNA_Node_arr_M1toM2 = []

    for node_name in unique_all:
        i = 0
        for seq in rnadata['mgi_symbol']:
            if node_name.lower() == seq.lower():
                dataM1 =   float(rnadata['M1'][i])
                dataM2 =   float(rnadata['M2'][i])
                dataM1M2 = float(rnadata['M1M2\n'][i])
                
                M2_mRNA_Node[node_name] = dataM2
                M2_mRNA_Node_arr.append(dataM2)
                M1_mRNA_Node[node_name] = dataM1
                M1_mRNA_Node_arr.append(dataM1)
                M1M2_mRNA_Node[node_name] = dataM1M2
                M1M2_mRNA_Node_arr.append(dataM1M2)
                
                # Convert any possible zero to small number
                ## division can screw this up 
                if dataM2 < 1e-14:
                    dataM2 = 1e-14

                if dataM1 < 1e-14:
                    dataM1 = 1e-14

                if dataM1M2 < 1e-14:
                    dataM1M2 = 1e-14

                # FC from M2 to M1M2
                logFCM12M2 = math.log2(dataM1M2/dataM2)
                logFC_mRNA_Node_M2toM12[node_name] = logFCM12M2 
                logFC_mRNA_Node_arr_M2toM12.append(logFCM12M2)
                # FC from M12 to M2
                logFC_mRNA_Node_M12toM2[node_name] = -logFCM12M2
                logFC_mRNA_Node_arr_M12toM2.append(-logFCM12M2)
                # FC from M2 to M1
                logFCM1M2 = math.log2(dataM1/dataM2)
                logFC_mRNA_Node_M2toM1[node_name] = logFCM1M2
                logFC_mRNA_Node_arr_M2toM1.append(logFCM1M2)
                # FC from M1 to M2
                logFC_mRNA_Node_M1toM2[node_name] = -logFCM1M2
                logFC_mRNA_Node_arr_M1toM2.append(-logFCM1M2)

            i += 1 

    ## Convert the mRNA based score on Node 
    M2_node_score = {}
    per10M2 = np.percentile(M2_mRNA_Node_arr,10)
    per50M2 = np.percentile(M2_mRNA_Node_arr,50)
    per95M2 = np.percentile(M2_mRNA_Node_arr,95)
    for key, value in M2_mRNA_Node.items():
        raw_value = float(value)
        if raw_value <= per10M2 :
            M2_node_score[key] = 1
        elif raw_value <= per50M2 and raw_value > per10M2:
            M2_node_score[key] = 2
        elif raw_value <= per95M2 and raw_value > per50M2:
            M2_node_score[key] = 3
        elif raw_value > per95M2:
            M2_node_score[key] = 4

    ## Convert the mRNA based score on Node 
    M1_node_score = {}
    per10M1 = np.percentile(M1_mRNA_Node_arr,10)
    per50M1 = np.percentile(M1_mRNA_Node_arr,50)
    per95M1 = np.percentile(M1_mRNA_Node_arr,95)
    for key, value in M1_mRNA_Node.items():
        raw_value = float(value)
        if raw_value <= per10M1 :
            M1_node_score[key] = 1
        elif raw_value <= per50M1 and raw_value > per10M1:
            M1_node_score[key] = 2
        elif raw_value <= per95M1 and raw_value > per50M1:
            M1_node_score[key] = 3
        elif raw_value > per95M1:
            M1_node_score[key] = 4

    ## Convert the mRNA based score on Node 
    M1M2_node_score = {}
    per10M1M2 = np.percentile(M1M2_mRNA_Node_arr,10)
    per50M1M2 = np.percentile(M1M2_mRNA_Node_arr,50)
    per95M1M2 = np.percentile(M1M2_mRNA_Node_arr,95)
    for key, value in M1M2_mRNA_Node.items():
        raw_value = float(value)
        if raw_value <= per10M1M2 :
            M1M2_node_score[key] = 1
        elif raw_value <= per50M1M2 and raw_value > per10M1M2:
            M1M2_node_score[key] = 2
        elif raw_value <= per95M1M2 and raw_value > per50M1M2:
            M1M2_node_score[key] = 3
        elif raw_value > per95M1M2:
            M1M2_node_score[key] = 4

    # Convert the mRNA based score on Node 
    # FC from M2 to M1M2
    FC_node_score_M2toM12 = {}
    for key, value in logFC_mRNA_Node_M2toM12.items():
        raw_value = float(value)
        if raw_value <= 0 :
            FC_node_score_M2toM12[key] = -1 # This is Bad
        else:
            FC_node_score_M2toM12[key] = 1 # This is Good

    # FC from M12 to M2
    FC_node_score_M12toM2 = {}
    for key, value in logFC_mRNA_Node_M12toM2.items():
        raw_value = float(value)
        if raw_value <= 0 :
            FC_node_score_M12toM2[key] = -1 # This is Bad
        else:
            FC_node_score_M12toM2[key] = 1 # This is Good

    # FC from M2 to M1
    FC_node_score_M2toM1 = {}
    for key, value in logFC_mRNA_Node_M2toM1.items():
        raw_value = float(value)
        if raw_value <= 0 :
            FC_node_score_M2toM1[key] = -1 # This is Bad
        else:
            FC_node_score_M2toM1[key] = 1 # This is Good

    # FC from M1 to M2
    FC_node_score_M1toM2 = {}
    for key, value in logFC_mRNA_Node_M1toM2.items():
        raw_value = float(value)
        if raw_value <= 0 :
            FC_node_score_M1toM2[key] = -1 # This is Bad
        else:
            FC_node_score_M1toM2[key] = 1 # This is Good

    # Scoring the edges #########################################
    '''
    This is where the magic happens
    ''' 

    edge_mRNA_score_M2toM12 = []
    edge_mRNA_score_M12toM2 = []
    edge_mRNA_score_M2toM1 = []
    edge_mRNA_score_M1toM2 = []
    
    for n in np.arange(num_of_edges):
        feature = networkdata['edge_features'][n]
        startnode = networkdata['start_nodes'][n]
        endnode = networkdata['end_nodes'][n]
        
        try: # by doing so, we are skipping some ligands that are not registered in mRNA seq **
            startM2  = M2_node_score[startnode]
            startM1  = M1_node_score[startnode]
            startM12 = M1M2_node_score[startnode]

        except: 
            startM2 = 0
            startM1 = 0 
            startM12 = 0 

        try: # by doing so, we are skipping some ligands that are not registered in mRNA seq **
            endM2    = M2_node_score[endnode]
            endM1    = M1_node_score[endnode]
            endM12   = M1M2_node_score[endnode]

        except: 
            #pass
            endM2 = 0
            endM1 = 0 
            endM12 = 0
        
        try: # by doing so, we are skipping some ligands that are not registered in mRNA seq **
            startFC_M2toM12 = FC_node_score_M2toM12[startnode]
            startFC_M12toM2 = FC_node_score_M12toM2[startnode]
            startFC_M2toM1  = FC_node_score_M2toM1[startnode]
            startFC_M1toM2  = FC_node_score_M1toM2[startnode]
            
        except: 
            startFC_M2toM12 = 0
            startFC_M12toM2 = 0
            startFC_M1toM2 = 0
            startFC_M2toM1 = 0 

        try: # by doing so, we are skipping some ligands that are not registered in mRNA seq **
            endFC_M2toM12   = FC_node_score_M2toM12[endnode]
            endFC_M12toM2   = FC_node_score_M12toM2[endnode]
            endFC_M2toM1    = FC_node_score_M2toM1[endnode]
            endFC_M1toM2    = FC_node_score_M1toM2[endnode]
            
        except: 
            endFC_M2toM12 = 0
            endFC_M12toM2 = 0
            endFC_M1toM2 = 0
            endFC_M2toM1 = 0

        if feature == 'up-regulates':
            score_M2toM12 = math.exp(-((startM2 + endM2)/2+startFC_M2toM12+endFC_M2toM12))
            score_M12toM2 = math.exp(-((startM12 + endM12)/2+startFC_M12toM2+endFC_M12toM2))
            score_M2toM1  = math.exp(-((startM2 + endM2)/2+startFC_M2toM1+endFC_M2toM1))
            score_M1toM2  = math.exp(-((startM1 + endM1)/2+startFC_M1toM2+endFC_M1toM2))
        elif feature == 'down-regulates':
            score_M2toM12 = math.exp(((startM2 + endM2)/2+startFC_M2toM12+endFC_M2toM12)) + 100 # I could add 1000 but let's see how it goes
            score_M12toM2 = math.exp(((startM12 + endM12)/2+startFC_M12toM2+endFC_M12toM2)) + 100 # I could add 1000 but let's see how it goes
            score_M2toM1  = math.exp(((startM2 + endM2)/2+startFC_M2toM1+endFC_M2toM1)) + 100 # I could add 1000 but let's see how it goes
            score_M1toM2  = math.exp(((startM1 + endM1)/2+startFC_M1toM2+endFC_M1toM2)) + 100 # I could add 1000 but let's see how it goes

        edge_mRNA_score_M2toM12.append(score_M2toM12)
        edge_mRNA_score_M12toM2.append(score_M12toM2)
        edge_mRNA_score_M2toM1.append(score_M2toM1)
        edge_mRNA_score_M1toM2.append(score_M1toM2)

    return unique_all, networkdata, edge_mRNA_score_M2toM12, edge_mRNA_score_M12toM2, edge_mRNA_score_M2toM1, edge_mRNA_score_M1toM2
```

## Generating network figure via NetworkX and perform analysis based on the weight assigned by the preceeding script
- The png file will be stored in your current working directory. 
- In this particular example, the current working directory is assigned to be /content/pathwayanalysis
- The file name is 'network_figure_XXXX.png'

```python
def networkAnalaysis(unique_all, networkdata, edge_mRNA_score, figure_file_name, receptorlistFilename, destination,listofNodes_path):
    ## Drawing part 
    G2 = Digraph('unix', filename='fullpicture',format='png',
            node_attr={'shape':'box',
                       'color': 'blue:purple', 
                       'style': 'filled',
                       'fontcolor':'white'},
            edge_attr={'color': 'red'})
  
    for node in unique_all:
      G2.node(node)

    for n in np.arange(len(networkdata['start_nodes'])):
      weight = edge_mRNA_score[n]
      G2.edge(networkdata['start_nodes'][n],networkdata['end_nodes'][n],label=str(weight))

    G2.view(figure_file_name)

    ## Actual Network Analysis
    G = nx.DiGraph()
    for node in unique_all:
      G.add_node(node)
    
    G.nodes()

    for n in np.arange(len(networkdata['start_nodes'])):
      weight = edge_mRNA_score[n]
      G.add_edge(networkdata['start_nodes'][n],networkdata['end_nodes'][n],weight=weight)   
    G.edges()

    Receptor_info = textreader(receptorlistFilename)
    Receptor_name = Receptor_info['Receptors']
    Receptor_feature = Receptor_info['Pro/Anti']
    Outcome = destination

    Scorelist = []
    ReceptorList = []
    ReceptorFeatureList = []

    PathwayNodes = {}

    i = 0
    for R in Receptor_name:
      try:
        Path = nx.dijkstra_path(G,R,Outcome)
        Score = nx.dijkstra_path_length(G,R,Outcome)
        Scorelist.append(Score)
        ReceptorList.append(R)
        ReceptorFeatureList.append(Receptor_feature[i])
        PathwayNodes[R] = Path
        if listofNodes_path:
          print('The shortest path from ',R,' (Upstream Receptor) to ',Outcome,' in a weighted Pathway Network: \n', Path,'\n with score of ',Score,'\n')
      except:
        pass
      i += 1

    collected_data = {'Receptors': ReceptorList, 'Pro/Anti': ReceptorFeatureList, 'Score': Scorelist}

    return collected_data, PathwayNodes
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
              ligand_list_file_name = 'ligandlist',
              destination = 'M1_polarization',
              unit_test_receptor = False,
              display=True,
              listofNodes_path=True,
              ): 
  ##
  ## Misc
  ## 
  
  unique_all, networkdata, edge_mRNA_score_M2toM12, edge_mRNA_score_M12toM2, edge_mRNA_score_M2toM1, edge_mRNA_score_M1toM2 = util.prepare_process(network_inputfilename,mRNA_inputfilename)
  
  if destination == 'M1_polarization':
    collected_data1, PathwayNodes1 = util.networkAnalaysis(unique_all, networkdata, edge_mRNA_score_M2toM12, 'new_figure_target'+destination, receptor_list, destination,listofNodes_path)
    collected_data2, PathwayNodes2 = util.networkAnalaysis(unique_all, networkdata, edge_mRNA_score_M2toM1, 'new_figure_control'+destination, receptor_list, destination,listofNodes_path)

  elif destination == 'M2_polarization':
    collected_data1, PathwayNodes1 = util.networkAnalaysis(unique_all, networkdata, edge_mRNA_score_M12toM2, 'new_figure_target'+destination, receptor_list, destination,listofNodes_path)
    collected_data2, PathwayNodes2 = util.networkAnalaysis(unique_all, networkdata, edge_mRNA_score_M1toM2, 'new_figure_control'+destination, receptor_list, destination,listofNodes_path)

  df_target  = pd.DataFrame(collected_data1)
  df_control = pd.DataFrame(collected_data2)

  network_inputfilename = re.sub('\.sif','',network_inputfilename)
  mRNA_inputfilename = re.sub('\.csv','',mRNA_inputfilename)

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

  with open(ligand_list_file_name+'.yaml') as file:
    ligand_list = yaml.load(file)

  receptor_specifics = {}
  ligand_specifics = {}

  for R in receptor_list['Receptors']:
    receptor_specifics[R] = {}
    ligand_specifics[R] = {}

    try:
      list_node = PathwayNodes1[R]
      for node in list_node:
        receptor_specifics[R][node] = {}
        prot_exp, M1, M2, M1M2, logFCofM1M2overM2 = util.expression_search(node,
                                                                           data_mrna,
                                                                           data_prot,
                                                                           prot_name)

        receptor_specifics[R][node]['mAb(From M1 to M2: M2/M1)'] = prot_exp
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

        ligand_specifics[R][ligand]['mAb(From M1 to M2: M2/M1)'] = prot_exp
        ligand_specifics[R][ligand]['M1_mRNA'] = M1
        ligand_specifics[R][ligand]['M2_mRNA'] = M2
        ligand_specifics[R][ligand]['M1M2_mRNA'] = M1M2
        ligand_specifics[R][ligand]['logFC(M1M2/M2)'] = logFCofM1M2overM2

    except:
      pass

  return collected_data1, collected_data2, df_target, df_control, receptor_specifics, ligand_specifics

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