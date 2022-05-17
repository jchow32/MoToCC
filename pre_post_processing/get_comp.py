# get_comp.py
# given input files of solution (cell \t value) and metadata of previous labels,
# join the files together and calculate and sort the percentage of each label


import sys
import pandas as pd


def main():
    sol_file = sys.argv[1]
    label_file = sys.argv[2]
    out_file = sys.argv[3]

    sol_df = pd.read_csv(sol_file, sep="\t")
    sol_df.columns = ['Cell', 'ObjectiveFunctionValue']
    label_df = pd.read_csv(label_file)[['Cell', 'Cluster']]
    # Supposing that there are columns called 'Cell' and 'Cluster' in the header of the label file
    cell_counts = pd.value_counts(label_df.Cluster).to_frame()

    df = pd.merge(sol_df, label_df, on='Cell', how='left').sort_values(by='Cluster')
    num_selected = df.shape[0]
    # print(num_selected)
    count_df = (pd.value_counts(df.Cluster)).to_frame()
    count_df['percent_composition'] = round((count_df.Cluster / num_selected) * 100, 2)
    # join with the cell_counts based on Cluster
    count_df = pd.merge(count_df.reset_index(), cell_counts.reset_index(), on='index', how='left')
    count_df['percent_of_cell'] = round((count_df.Cluster_x / count_df.Cluster_y) * 100, 2)
    # count_df = count_df.set_index('index')
    count_df.columns = ['cell-type', 'num_selected', 'per_comp', 'total_cell', 'per_capture']
    print(count_df)
    count_df.to_csv(out_file, sep="\t", mode='a', header=False)
    # user will have to add their own header to this file after all appending is finished


if __name__ == "__main__":
    main()
