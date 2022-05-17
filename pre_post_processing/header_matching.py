# header_matching.py
# given the header of a file (comma-delimited),
# and a key of corresponding geneID \t geneSymbol
# substitute the header with the corresponding geneSymbol and write the new header to file
# use bash to cat together the tail -n+2 of the exp file with new header.


import sys
import pandas as pd


def main():
	exp_file = sys.argv[1]
	gene_key = sys.argv[2]
	out_file = sys.argv[3]

	exp_df = pd.read_csv(exp_file, nrows=1)
	key_df = pd.read_csv(gene_key, sep="\t")

	# turn key_df into dictionary, where the key is the stable ID and the value is the gene symbol
	# first, keep only the keys which DON'T have NaN for gene symbol.
	key_df = key_df[key_df['Gene name'].notna()]
	key_dict = dict(key_df.values)

	# now replace the corresponding columns in exp_df with those from key_dict
	exp_df = exp_df.rename(columns=key_dict)
	exp_df.to_csv(out_file, sep=",", index=False)


if __name__ == "__main__":
	main()
