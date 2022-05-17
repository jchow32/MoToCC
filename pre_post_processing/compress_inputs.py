# compress_inputs.py
# compress the expression, KNN, and SNN files using hdf5


import sys
import pandas as pd
from scipy import sparse


def wide_df_to_hdf(filename, data, columns=None, maxColSize=2000, **kwargs):
    # https://stackoverflow.com/questions/16639503/unable-to-save-dataframe-to-hdf5-object-header-message-is-too-large
    """Write a `pandas.DataFrame` with a large number of columns
    to one HDFStore.

    Parameters
    -----------
    filename : str
        name of the HDFStore
    data : pandas.DataFrame
        data to save in the HDFStore
    columns: list
        a list of columns for storing. If set to `None`, all
        columns are saved.
    maxColSize : int (default=2000)
        this number defines the maximum possible column size of
        a table in the HDFStore.

    """
    import numpy as np
    from collections import ChainMap
    store = pd.HDFStore(filename, **kwargs)
    if columns is None:
        columns = data.columns
    colSize = columns.shape[0]
    if colSize > maxColSize:
        numOfSplits = np.ceil(colSize / maxColSize).astype(int)
        colsSplit = [
            columns[i * maxColSize:(i + 1) * maxColSize]
            for i in range(numOfSplits)
        ]
        _colsTabNum = ChainMap(*[
            dict(zip(columns, ['data{}'.format(num)] * colSize))
            for num, columns in enumerate(colsSplit)
        ])
        colsTabNum = pd.Series(dict(_colsTabNum)).sort_index()
        for num, cols in enumerate(colsSplit):
            store.put('data{}'.format(num), data[cols], format='table')
        store.put('colsTabNum', colsTabNum, format='fixed')
    else:
        store.put('data', data[columns], format='table')
    store.close()


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
            # the index is where it comes from, and the data is the gene name?
            # colsTabNum.index.tolist() ??
            # tabNums = pd.Series(index=colsTabNum[columns].values, data=colsTabNum[columns].data).sort_index()
            tabNums = pd.Series(data=colsTabNum[columns].values, index=colsTabNum[columns].index.tolist()).sort_index()

            # i think it's the other way around?
            for table in tabNums.unique():
                # table is the unique dataX, tabNums[table] is the list of genes
                # growing_data.append(store.select(table, columns=tabNums[table], **kwargs))
                growing_data.append(store.select(table, columns=tabNums.index.tolist()))
        else:
            for table in colsTabNum.unique():
                growing_data.append(store.select(table, **kwargs))
        growing_data = pd.concat(growing_data, axis=1).sort_index(axis=1)
    else:
        growing_data = store.select('data', columns=columns)
    store.close()
    return growing_data


def write_sparse(in_file, out_string, out_name):
    print("Reading input file")
    df = pd.read_csv(in_file)
    print("Writing sparse output")
    out_parse = sparse.csr_matrix(df.values)
    sparse.save_npz("%s/%s" % (out_string, out_name), out_parse)


def main():
    exp_file = sys.argv[1]
    sim_file = sys.argv[2]
    knn_file = sys.argv[3]
    out_string = sys.argv[4]

    print("Reading expression file")
    exp_df = pd.read_csv(exp_file)
    wide_df_to_hdf("%s/exp_norm_scale.hdf5" % out_string, exp_df)
    # df = read_hdf_wide_df("%s/exp_norm_scale.hdf5" % out_string)

    write_sparse(sim_file, out_string, "sim_sparse")
    write_sparse(knn_file, out_string, "knn_sparse")


if __name__ == "__main__":
    main()
