# Cell-type resolved transcriptional network analysis of in vivo cellular senescence following injury

This repository contains the code corresponding to the paper "Cell-type resolved transcriptional network analysis of in vivo cellular senescence following injury", by Alda Sabalic, Victoria Moiseeva, Andres Cisneros, Oleg Deryagin, Eusebio Perdiguero, Pura Muñoz-Cánoves, and Jordi Garcia-Ojalvo.

The main pipeline consists of a collection of Python notebooks and scripts, containing functions to filter the original dataset, compute principal component analysis (PCA) of the filtered data, and combine an eigenvector centrality analysis of a gene discrimination network with community detection of a cell-condition network. For details of the method, see the paper above. The eigenvector centrality analysis was performed in Matlab, using code developed by Roffo and Melzi ("Ranking to learn: Feature ranking and selection via eigenvector centrality". In International Workshop on New Frontiers in Mining Complex Patterns 2016, pp. 19-35, Springer).

The repository also includes a validation using single-cell RNA sequencing data published previously by Omori et al ("Generation of a p16 reporter mouse and its use to characterize and target p16high cells in vivo". Cell Metabolism 32, pp. 814-828, 2020).
