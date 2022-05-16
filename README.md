# C-TULPSE
Cell-types using linear programming and selective expression (C-TULPSE) is a linear programming approach to identify cells that selectively express a set of genes of interest, such as the genes belonging to a genetic module. Given single-cell or single-nucleus gene expression data, cell-cell similarity and K-nearest neighbor connectivity derived from the expression data, module genes, and an upper bound to the number of cells to return as a solution (<i>k</i>), C-TULPSE first calculates edge weights among module genes as per Equation 1, then uses linear programming to find an initial solution consisting of approximately <i>k</i> cells. The initial solution of candidate cells is refined by identifying cells in the largest strongly connected component (LSCC) and returning the set of candidate cells in the LSCC as the solution (Figure 1). 

INSERT_FIGURE_1

C-TULPSE maximizes an objective function (Equation 1) defined as the product of the correlated co-expression of module genes and cell-cell similarity while subject to linear constraints (Equation 2). 

INSERT_OBJECTIVE_FUNCTION 
INSERT_LINEAR_CONSTRAINTS

# How to run C-TULPSE
Python packages needed to run C-TULPSE are listed in LINK_TO_CONDA_PACKAGE. <br>

  Required parameters: <br>
  - $expression <normalized expression data (.hdf5)>
    - The expression file contains normalized and scaled single-cell expression data (cells x genes)
    - LINK_TO_SCRIPT compresses expression data 
  - $similarity <cell-cell similarity derived from expression data (.npz)>
    - Cell-cell similarity (cells x cells) can be retrieved as the shared-nearest neighbor graph via the script LINK_TO_R_SCRIPT
    - LINK_TO_SCRIPT saves cell-cell similarity as a sparse matrix
  - $knn <K-nearest neighbor connectivity derived from expression data (.npz)>
    - K-nearest neighbor (KNN) connectivity can be retrieved via the script LINK_TO_R_SCRIPT 
    - LINK_TO_SCRIPT saves the KNN graph as a sparse matrix (cells x cells)
  - $module <genes of interest>
    - Text file consisting of module genes (new-line delimited)
  - $cell_ids <cell ids from expression data>
    - Text file containing all cell ids in the **same** order as in $expression (new-line delimited)
    - LINK_TO_SCRIPT retrieves cell_ids from expression data
  - $gene_ids <gene ids from expression data>
    - Text file containing all gene ids in the **same** order as in $expression (new-line delimited)
    - LINK_TO_SCRIPT retrieves gene_ids from expression data
  - $k <upper bound of cells to return as solution>
    - Integer
  - $output_string <output string to contain in output file names>
    - String
  - $quickstart <string to enable or disable quickstart>
    - String to enable (True) or disable (False) quickstart
    - Quickstart should be disabled the first time a module is provided to C-TULPSE (see Tips)
    
Example command with quickstart disabled: 
```
python3 c_tulpse.py \
$expression  \
$similarity  \
$knn \
$module \
$cell_ids \
$gene_ids \
$k \
$output_string \
False
```

# Tips
* The parameter _k_ can be varied for any given module to select cells at varying degrees of resolution. A silhouette score describing the similiarty of selected cells following two-dimensional t-SNE reduction can be calculated via LINK_TO_SCRIPT.
* Percent composition of selected cells can be calculated via LINK_TO_SCRIPT if cell-type labels are known.
* After running C-TULPSE for a given module for the first time (quickstart=False), edge weights associated with the module are automatically saved. If the user wishes to provide multiple, unique values of _k_ to C-TULPSE, enabling the quickstart option allows C-TULPSE to read the saved edge weights and reduce runtime on subsequent runs.
* Example pre- and post-processing scripts are available for the compression of expression data, percent composition calculation, and visualization of selected cells via dimensionality reduction plots. 
