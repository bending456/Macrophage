---
title: Cytoscape
nav_order: 1
parent: User Guides
---

# User Guides

## How to Start
1. Download SIF file from the repository: [Click to Download](https://bitbucket.org/pkhlab/pathwayanalysis/get/e2b9464490d5.zip) 
2. Unzip the download file 
3. Locate .sif file (*expanded_IL6_modified_receptor.sif* for this tutorial)

### In Cytoscape 
- **[Main Tab]** File -> Import -> Network from File or [Ctrl + L]
![This is the first visual in cytoscape](cyto1.png)
- **[Main Tab]** Layout -> yFiles Hierarchic Layout (This is a plugin so you may be required to download and install as yFiles package)
    - To install yFiles Layout Algorithms: **[Main Tab]** Apps -> App Manager -> Search yfiles in the searching bar -> Install 
    ![This is app manager](cyto2.png)

![This is refined network in cytoscape](cyto3.png)

#### Exploring the network 
- **[Search Bar]** Type the target node (example: IL1R1) and press [Enter]

![This is search bar](cyto4.png)

- The highlighted node and edges are associated with searched node
![This is searched outcome](cyto5.png)

- **[Zoomed with scroll]**

![This is zoomed outcome](cyto6.png)

### Customizing Visualization 
- **[Side Tab]** Style -> Edge -> Middle block of [*Target Arrow Shape*] -> Discrete Mapping in [*Mapping Type*]

![This is side tab](cyto7.png)

- **[Zoomed with Scroll]**

![This is zoomed](cyto8.png)