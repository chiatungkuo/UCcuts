## Unified and Constrasting Cuts (UCcuts) in Multiple Graphs

This repository contains source code to learn a unified graph cut from a collection of graphs and a constrasting cut between two collections of graphs. Graphs in the collections are assumed to be of the same modality and in some sense "compatible".
The details of the methods are described in my paper [Unified and Contrasting Cuts in Multiple Graphs: Application to Medical Imaging Segmentation](https://sites.google.com/site/chiatungkuo/publication), which demonstrates their utility in the analysis of brain scans. The data used in the paper was from UC Davis Alzheimer's Institute; its access may be granted upon request (and appropriate disclosure agreement). Public data of similar format could also be obtained from the [Human Connectome Project](http://www.humanconnectome.org/data/).

The source code also includes my implementation of one popular multivew spectral clustering method (Zhou et. al. 2012) and a consensus spectral clustering method adapted from another recent work (Lancichinetti and Fortunato, 2012). 

File(s) in this repository: 

+ unifiedcut.m: learn a unified cut from a set of graphs
+ contrastcut.m: learn a contrasting cut between two sets of graphs
+ consensus.m: learn a consensus cut from a set of graphs (adapted from work by Lancichinetti and Fortunato, 2012)
+ multisp.m: multiview spetral clustering method (implementation of work by Zhou et. al., 2012)

Please report any bugs/issues to [tomkuo@ucdavis.edu](mailto:tomkuo@ucdavis.edu).
 
