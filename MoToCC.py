# MoToCC.py
# given normalized expression information for cells and their similarity, choose an optimal subset that maximizes
# the objective function

# given a graph consisting of cells (c) in C and edges (i, j) in E,
# Where z_i,j is the weighted pairwise autocorrelation between two cells i and j,
# x_i,j describes whether an edge is selected or not between cells i and j,
# k is the number of cells to select,
# genes (g) from gene n to p in a module G,
# and w_i,j is the cell-cell similarity between cells i and j,
# and n_i is the normalized expression of gene n for cell i, etc.

# maximize the summation of z_i,j * x_i,j
# where z_i,j is the summation (summation) for all genes in G, w_i,j * (n_i*p_j + p_i*n_j)

# with constraints
# 0 <= x_i,j <= 1 for all (i,j) in E -- ADDED
# summation y_i <= k  -- ADDED
# x_i,j <= y_j -- ADDED
# x_i,j <= y_i -- ADDED
# x_i,j >= y_i + y_j - 1  -- ADDED

import sys
from docplex.mp.model import Model
import pandas as pd
import numpy as np
import time
import itertools
import networkx as nx
from scipy import sparse
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import argparse


def get_non_zero_sim_cells(sim_df, cells_list):
    idx_values = sparse.find(sim_df)
    rows_idx = idx_values[0]
    cols_idx = idx_values[1]
    non_zero_sim_idx = np.vstack((rows_idx, cols_idx)).transpose()

    non_zero_sim_idx = [subl for subl in non_zero_sim_idx if subl[1] != subl[0]]
    # remove repetitive combinations, like [1, 1], [2,2]
    non_zero_sim_idx = set(map(lambda x: tuple(sorted(x)), non_zero_sim_idx))
    # remove repetitive combinations like [1,2] vs [2,1]

    # create a dictionary which has the cell as the key and the index (start at 0) as the value
    cell_dict = {k: v for k, v in enumerate(cells_list)}
    # cells that have non-zero similarity between them- write a list that replaces the indices
    # with their corresponding cell labels
    valid_cell_list = ([(cell_dict.get(cell1), cell_dict.get(cell2)) for cell1, cell2 in non_zero_sim_idx])
    # print(valid_cell_list)  # it's a list of tuples of pairs of cells that have non-zero sim. with each other

    return valid_cell_list, non_zero_sim_idx, cell_dict


def read_hdf_wide_df(filename, columns=None, **kwargs):
    # https://stackoverflow.com/questions/16639503/unable-to-save-dataframe-to-hdf5-object-header-message-is-too-large
    """Read a `pandas.DataFrame` from a HDFStore.

    Parameter
    ---------
    filename : str
        name of the HDFStore
    columns : list
        the columns in this list are loaded. Load all columns,
        if set to `None`.

    Returns
    -------
    data : pandas.DataFrame
        loaded data.

    """
    store = pd.HDFStore(filename)
    growing_data = []
    colsTabNum = store.select('colsTabNum')
    if colsTabNum is not None:
        if columns is not None:
            tabNums = pd.Series(data=colsTabNum[columns].values, index=colsTabNum[columns].index.tolist()).sort_index()
            for table in tabNums.unique():
                growing_data.append(store.select(table, columns=tabNums.index.tolist()))
        else:
            for table in colsTabNum.unique():
                growing_data.append(store.select(table, **kwargs))
        growing_data = pd.concat(growing_data, axis=1).sort_index(axis=1)
    else:
        growing_data = store.select('data', columns=columns)
    store.close()
    return growing_data


def read_files(exp_file, input_module, sim_file, knn_file, cell_file, gene_file):
    module_file_name = input_module.split("/")[-1]
    weight_file = '/'.join(exp_file.split('/')[0:-1]) + ("/edge_weights_%s" % module_file_name)
    # note that the edge_weights are EXACTLY the same regardless of k

    with open(input_module) as module_file:
        module_genes = module_file.read().splitlines()  # store the module genes into a list
    module_combos = list(itertools.combinations(range(len(module_genes)), 2))  # get all non-redundant pairs of genes

    with open(cell_file) as cells:
        cells_list = cells.read().splitlines()  # store the module genes into a list

    with open(gene_file) as genes:
        genes_list = genes.read().splitlines()  # store the module genes into a list

    # print out what genes are in module_genes but not avail_genes for a re-run
    missing_genes = list(set(module_genes) - set(genes_list))
    if len(missing_genes) >= 1:
        print("Genes in provided module that do not exist in dataset:")
        for gene in missing_genes:
            print(gene)
        print("Remove or rename listed genes from module. Exiting")
        exit()
    module_genes = list(set(module_genes).intersection(genes_list))

    exp_df = read_hdf_wide_df(exp_file, columns=module_genes)

    sim_df = sparse.load_npz(sim_file)
    knn_df = sparse.load_npz(knn_file)

    return weight_file, module_genes, module_combos, exp_df, sim_df, knn_df, cells_list, genes_list


def prune_edges(weight_file, valid_cell_list, prune):
    # setting single scores equal to 0 when they within 1 SD of the mean score
    weight_df = pd.read_csv(weight_file, header=None, sep="\t")  # read the file you just wrote of the edge weights
    weight_df.columns = ['cell_pair', 'edge_weight']
    # of the format cell_cell \t weight
    # it's a list of tuples of pairs of cells that have non-zero sim. with each other

    w_mean = weight_df['edge_weight'].mean()
    w_std = weight_df['edge_weight'].std()
    one_up = w_mean + w_std
    one_down = w_mean - w_std

    if prune:
        keep = (weight_df[(weight_df['edge_weight'] >= one_up) | (weight_df['edge_weight'] <= one_down)])
        keep[['cell_1', 'cell_2']] = keep['cell_pair'].str.split('_', expand=True)
        valid_cell_list = list(zip(keep.cell_1, keep.cell_2))
    # retain the non-zero score for those outside one SD

    weight_dict = dict(weight_df.values)  # convert the data frame into a dictionary for maximization
    # although the weight_dict contains entries that do not exist in modified valid_cell_list, they won't be
    # referenced anyways

    return weight_dict, valid_cell_list


def calculate_gene_scores(sim_df, exp_df, valid_cell_idx, module_combos, weight_file, cell_dict):

    exp_array = exp_df.as_matrix()

    for cell_pair in valid_cell_idx:
        edge_weight = 0
        i = cell_pair[0]  # this is the label of a single cell, connected to j
        j = cell_pair[1]  # this is the label of a single cell, connected to i
        weight = sim_df[i, j]

        # need non-redundant pairs of module genes, excluding connection with self
        for gene_pair in module_combos:
            gene_1 = gene_pair[0]
            gene_2 = gene_pair[1]
            n_i = exp_array[i][gene_1]
            p_j = exp_array[j][gene_2]
            p_i = exp_array[i][gene_1]
            n_j = exp_array[j][gene_2]

            edge_weight += weight * ((n_i * p_j) + (p_i * n_j))

        # when you write the edges, write i and j as their cell names instead of indices. Use the dictionary for this
        with open(weight_file, 'a') as edge_file:
            edge_file.write("%s_%s\t%s\n" % (cell_dict[i], cell_dict[j], edge_weight))


def get_strong_cc(knn_df, selected_cells, cell_dict):

    reverse_cell_dict = {v: k for k, v in cell_dict.items()}
    subset_cell_idx_list = sorted(list({your_key: reverse_cell_dict[your_key] for your_key in selected_cells}.values()))
    # this is a list of the indices of the selected cells

    # you can directly slice the sparse matrix because you know the indices of the selected cells
    knn_subset = knn_df[subset_cell_idx_list, :][:, subset_cell_idx_list]

    # you must enumerate the subset_cell_idx_list
    k_dict = {k: v for k, v in enumerate(subset_cell_idx_list)}
    # the key is 0, 1, 2, etc. with the value being the subset_cell_idx_list

    g = nx.from_scipy_sparse_matrix(knn_subset, create_using=nx.DiGraph())
    comp_list = {}
    largest_comp = 0

    for component in nx.strongly_connected_components(g):  # instead of connected_components
        # take the largest STRONGLY connected comp
        if len(component) > largest_comp:
            largest_comp = len(component)
            comp_list = list(component)
    # get the IDs of the cells that are in the biggest comp

    # in comp_list are the indices within the KNN subset. So you need to translate them into full data set ids
    certain_cells = [k_dict.get(key) for key in comp_list]  # returns indices you search for in cell_dict
    # get the keys that match certain_cells from the keys of cell_dict
    certain_cells = [cell_dict[key] for key in certain_cells]

    return certain_cells


def parse_solution(solution):
    r_solution_list = []
    solution_list = []
    for var, value in solution.iter_var_values():
        var = str(var).replace('c_', '').replace(")", '').replace("(", '').replace("'", '').replace(", ", "_")
        if "_" not in var:
            r_solution_list.append(var)
            solution_list.append([var, value])
    solution_df = pd.DataFrame(data=solution_list, columns=['cell', 'value'])  # keep track of obj values

    return solution_df


def linear_program(valid_cell_list, cells_list, num_cells, weight_dict, out_file):

    m = Model(name='critical_cell', checker='off')  # initialize cplex model

    indicator = m.continuous_var_dict(valid_cell_list, name='c', lb=0, ub=1)
    cell_indicator = m.continuous_var_dict(cells_list, name='c', lb=0, ub=1)

    m.add_constraints(indicator[cell_pair] <= cell_indicator[cell_pair[0]] for cell_pair in valid_cell_list)
    m.add_constraints(indicator[cell_pair] <= cell_indicator[cell_pair[1]] for cell_pair in valid_cell_list)
    m.add_constraints(indicator[cell_pair] >=
                      (cell_indicator[cell_pair[0]] + cell_indicator[cell_pair[1]] - 1)
                      for cell_pair in valid_cell_list)
    num_picked_cells = 0
    for cell in cells_list:
        num_picked_cells += cell_indicator[cell]
    m.add_constraint(num_picked_cells <= num_cells)  # only select the user-specified number of cells
    m.maximize(m.sum(
        weight_dict["%s_%s" % (chosen_cell[0], chosen_cell[1])] * indicator[
            chosen_cell[0], chosen_cell[1]] for chosen_cell in valid_cell_list))
    m.print_information()
    solution = m.solve()
    solution_df = parse_solution(solution)
    solution_df.to_csv(out_file, sep="\t", index=False)

    obj_value = m.objective_value
    print("Obj. function value: %s" % obj_value)


def main():
    parser = argparse.ArgumentParser()
    req_grp = parser.add_argument_group(title='Required')
    req_grp.add_argument("--exp", "-e", help="Normalized, scaled gene expression (.hdf5)", required=True)
    req_grp.add_argument("--sim", "-sim", help="Cell-cell similarity (.npz)", required=True)
    req_grp.add_argument("--knn", "-knn", help="K-nearest neighbor (.npz)", required=True)
    req_grp.add_argument("--module_genes", "-m", help="Module genes", required=True)
    req_grp.add_argument("--cell_ids", "-c", help="Cell IDs", required=True)
    req_grp.add_argument("--gene_ids", "-g", help="Gene IDs", required=True)
    req_grp.add_argument("--k", "-k", help="Upper bound cells in solution", required=True)
    req_grp.add_argument("--out_string", "-o", help="Output string", required=True)
    parser.add_argument(
        "--quickstart", help="Enable quickstart if edge weights previously calculated", action='store_true')
    parser.add_argument(
        "--prune", help="Enable pruning of edge weights", action='store_true')
    parser.add_argument(
        "--integer", help="Enable integer programming rather than linear programming", action='store_true')

    args = parser.parse_args()
    exp_file = args.exp
    sim_file = args.sim
    knn_file = args.knn
    input_module = args.module_genes
    cells = args.cell_ids
    genes = args.gene_ids
    num_cells = int(args.k)
    custom_string = args.out_string

    quick_start = args.quickstart
    prune = args.prune
    integer_bool = args.integer

    start_time = time.time()
    print("Reading input files: %s" % (time.time() - start_time))
    weight_file, module_genes, module_combos, exp_df, sim_df, knn_df, cells_list, genes_list = read_files(
        exp_file, input_module, sim_file, knn_file, cells, genes)
    print("Read input files: %s" % (time.time() - start_time))

    # reduce the number of variables in model by considering pairs of cells that have similarity > 0
    valid_cell_list, valid_cell_idx, cell_dict = get_non_zero_sim_cells(sim_df, cells_list)
    print("Returned pairs of cells with non-zero similarity: %s" % (time.time() - start_time))

    if not quick_start:  # run the model from scratch and calculate edge weights for first time
        print("Calculating gene scores: %s" % (time.time() - start_time))
        calculate_gene_scores(sim_df, exp_df, valid_cell_idx, module_combos, weight_file, cell_dict)
        del exp_df
        del sim_df
        print("Finished gene scoring calculations: %s" % (time.time() - start_time))

    weight_dict, valid_cell_list = prune_edges(weight_file, valid_cell_list, prune)

    m = Model(name='critical_cell', checker='off')  # initialize cplex model

    indicator = m.continuous_var_dict(valid_cell_list, name='c', lb=0, ub=1)
    cell_indicator = m.continuous_var_dict(cells_list, name='c', lb=0, ub=1)

    print("Loading batch constraints: %s" % (time.time() - start_time))
    m.add_constraints(indicator[cell_pair] <= cell_indicator[cell_pair[0]] for cell_pair in valid_cell_list)
    m.add_constraints(indicator[cell_pair] <= cell_indicator[cell_pair[1]] for cell_pair in valid_cell_list)
    m.add_constraints(indicator[cell_pair] >=
                      (cell_indicator[cell_pair[0]] + cell_indicator[cell_pair[1]] - 1)
                      for cell_pair in valid_cell_list)
    num_picked_cells = 0
    for cell in cells_list:
        num_picked_cells += cell_indicator[cell]
    m.add_constraint(num_picked_cells <= num_cells)  # only select the user-specified number of cells
    print("Loaded constraints: %s" % (time.time() - start_time))
    print("Beginning maximization of objective function: %s" % (time.time() - start_time))
    m.maximize(m.sum(
        weight_dict["%s_%s" % (chosen_cell[0], chosen_cell[1])] * indicator[
            chosen_cell[0], chosen_cell[1]] for chosen_cell in valid_cell_list))
    m.print_information()
    solution = m.solve()
    solution_df = parse_solution(solution)

    obj_value = m.objective_value
    print("Initial obj. function value: %s" % obj_value)
    
    s_file = '/'.join(exp_file.split('/')[0:-1]) + "/%s_%s_solution.txt" % (custom_string, num_cells)
    i_file = '/'.join(exp_file.split('/')[0:-1]) + "/%s_%s_init_solution.txt" % (custom_string, num_cells)
    sl_file = '/'.join(exp_file.split('/')[0:-1]) + "/%s_%s_LPonLSCC.txt" % (custom_string, num_cells)

    # this initial solution will consist of less than num_cells, and only cells with values > 0
    if integer_bool:
        certain_df = solution_df.sort_values(
            by='value', ascending=False).head(num_cells)
        certain_df.to_csv(s_file, index=False, sep="\t")
    else:
        # save the initial solution for double checking
        certain_df = solution_df.sort_values(
            by='value', ascending=False)
        certain_df.to_csv(i_file, index=False, sep="\t")

        selected_cells = certain_df['cell'].tolist()

        # find LSCC. Provide cells in the LSCC to the LP with the same constraints
        certain_cells = get_strong_cc(knn_df, selected_cells, cell_dict)
        # filter valid_cell_list to only contain cells in certain_cells
        new_valid_cell_list = []
        for valid_pair in valid_cell_list:
            cell_1 = valid_pair[0]
            cell_2 = valid_pair[1]
            if cell_1 in certain_cells and cell_2 in certain_cells:
                new_valid_cell_list.append(valid_pair)
        linear_program(new_valid_cell_list, certain_cells, num_cells, weight_dict, sl_file)

        weight_df = pd.read_csv(weight_file, header=None, sep="\t")  # read the file you just wrote of the edge weights
        weight_df.columns = ['cell_pair', 'edge_weight']

        # you need to filter weight_df so that only cells in certain_cells are listed
        filtered_edges = []
        for c_pair in weight_dict:
            split_c = c_pair.split('_')
            if split_c[0] in certain_cells and split_c[1] in certain_cells:
                filtered_edges.append(c_pair)

        weight_df = weight_df[weight_df['cell_pair'].isin(filtered_edges)]
        weight_df[['cell_1', 'cell_2']] = weight_df['cell_pair'].str.split('_', expand=True)
        weight_pivot = weight_df.pivot(index='cell_1', columns='cell_2', values='edge_weight').fillna(0)

        weight_pivot_embedded = TSNE().fit_transform(weight_pivot)
        kmeans = KMeans(n_clusters=2).fit_predict(weight_pivot_embedded)
        score = silhouette_score(weight_pivot_embedded, kmeans)
        print("Silhouette score: %s" % score)

    print("Time elapsed: %s" % (time.time() - start_time))


if __name__ == "__main__":
    main()
