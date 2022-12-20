# heuristic_summation.py
# Given individual cells and their labels, calculate s_h
# where s_h is the summation of (all cells i in C_h) summation (all genes in M) [g_i] / | C_h |, where
# C_h is the number of cells in cell-type h

import pandas as pd
import sys
import numpy as np


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


def read_module_calculate_sh(module_file, exp_file, meta_df, meta_dict, out_string):
    with open(module_file) as module_in:
        module_genes = module_in.read().splitlines()  # store the module genes into a list
    exp_df = read_hdf_wide_df(exp_file, columns=module_genes)

    sh_dict = {}

    for cell_type in meta_dict:
        # get the list of cells that are in this cell_type
        cells_list = meta_df[meta_df['Cluster'] == cell_type]['Cell'].tolist()  # a list of the cells for this cell-type
        # filter the exp_df for cells in this subset
        exp_subset_df = exp_df[exp_df.index.isin(cells_list)]
        # divide the matrix by the number of cells in this cell-type
        num_cells = meta_dict[cell_type]
        exp_subset_df = exp_subset_df / num_cells
        # sum over the rows, then the columns
        s_h = exp_subset_df.sum(axis=1).sum(axis=0)

        sh_dict[cell_type] = abs(s_h)

    print(sh_dict)
    most_rel = max(sh_dict, key=sh_dict.get)  # the key with max value
    biggest_sh = sh_dict[most_rel]
    print("%s: %s = %s" % (out_string, most_rel, biggest_sh))


def calculate_sh(mod_string, metadata_file, module_file, exp_file):
    # a list of the random module files
    rand_list = ["rand_%s_module_%s.txt" % (mod_string, x) for x in range(1, 21)]

    # for every cell-type in the dataset, for that cell-type determine HOW MANY cells exist, then
    # gather the normalized gene expression values for every gene in module for that cell i
    meta_df = pd.read_csv(metadata_file)[['Cell', 'Cluster']]  # comma sep
    # store the cell-type and the number of each type in a dictionary
    meta_dict = meta_df.groupby('Cluster')['Cell'].nunique().to_dict()

    read_module_calculate_sh(module_file, exp_file, meta_df, meta_dict, mod_string)
    # for i in rand_list:
    #     read_module_calculate_sh(i, exp_file, meta_df, meta_dict, i)
    # also calculate the average s_h for 20 randomized modules of the same size
    # they follow the format of rand_${mod}_module_${i}.txt


def calculate_sh_certain(sel_cells, exp_df):
    cell_list = pd.read_csv(sel_cells, sep="\t")['cell'].tolist()
    exp_subset_df = exp_df[exp_df.index.isin(cell_list)]
    num_cells = len(cell_list)
    exp_subset_df /= num_cells
    s_h = exp_subset_df.sum(axis=1).sum(axis=0)
    return (s_h)


def calculate_sh_k(exp_df, sel_cells):
    k = len(pd.read_csv(sel_cells, sep="\t")['cell'].tolist())  # how many cells selected by solver
    new_df = exp_df / k
    sel_cells = new_df.sum(axis=1).sort_values(ascending=False).head(k).index.tolist()

    exp_subset_df = exp_df[exp_df.index.isin(sel_cells)]
    exp_subset_df /= k
    s_h = exp_subset_df.sum(axis=1).sum(axis=0)
    return s_h


def main():
    metadata_file = sys.argv[1]
    module_file = sys.argv[2]
    exp_file = sys.argv[3]
    mod_string = sys.argv[4]

    # calculate_sh(mod_string, metadata_file, module_file, exp_file)
    with open(module_file) as module_in:
        module_genes = module_in.read().splitlines()  # store the module genes into a list
    exp_df = read_hdf_wide_df(exp_file, columns=module_genes)

    sh_dict = {}
    naive_k_dict = {}
    for k in range(250, 5250, 250):
        sel_cells = "%s_%s_LPonLSCC.txt" % (mod_string, k)
        # s_h = calculate_sh_certain(sel_cells, exp_df)
        # sh_dict[k] = s_h

        s_h = calculate_sh_k(exp_df, sel_cells)

        naive_k_dict[k] = s_h
    naive_out = pd.DataFrame(list(naive_k_dict.items()))
    naive_out.to_csv("%s_sh_naive_k.tsv" % mod_string, sep="\t", index=False)

    # sh_out = pd.DataFrame(list(sh_dict.items()))
    # sh_out.to_csv("%s_sh_certain_abs.tsv" % mod_string, sep="\t", index=False)

    # write to file
    # do similar for random
    """
    rand_sh_dict = {}
    sh_list = []
    for k in range(250, 5250, 250):
        print(k)
        for n in range(1, 21):
            # try if exist
            sel_cells = "Rand_%s_%s_%s_LPonLSCC.txt" % (mod_string, n, k)
            try:
                cell_list = pd.read_csv(sel_cells, sep="\t")['cell'].tolist()
            except FileNotFoundError:
                break

            rand_module_file = "rand_%s_module_%s.txt" % (mod_string, n)
            with open(rand_module_file) as module_in:
                module_genes = module_in.read().splitlines()  # store the module genes into a list
            exp_df = read_hdf_wide_df(exp_file, columns=module_genes)

            exp_subset_df = exp_df[exp_df.index.isin(cell_list)]

            num_cells = len(cell_list)
            exp_subset_df /= num_cells
            s_h = exp_subset_df.sum(axis=1).sum(axis=0)

            sh_list.append(s_h)
        if len(sh_list) > 0:
            avg_sh_for_k = np.mean(sh_list)
            rand_sh_dict[k] = avg_sh_for_k
        sh_list = []
    sh_rand_out = pd.DataFrame(list(rand_sh_dict.items()))
    sh_rand_out.to_csv("%s_sh_certain_rand.tsv" % mod_string, sep="\t", index=False)
    """


if __name__ == "__main__":
    main()
