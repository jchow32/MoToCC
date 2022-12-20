# proportion_gene_exp.py
# given a module and a starter string of the genes ranked by expression,
# find the proportion of each gene rank that exists in the true module,
# and then return a random list of the same size of true module, with correct proportions

import sys
import random


def main():
    module_file = sys.argv[1]
    gene_rank_string = sys.argv[2]  # cortex_top_exp_genes
    out_string = sys.argv[3]

    with open(module_file) as m_file:
        module_genes = m_file.read().splitlines()  # store the module genes into a list

    # search for intersection of module_genes with gene_rank_string files 1-10, and print out the shape of this
    # Do this 20 times for 20 randomized modules of same size as true
    for j in range(1, 21):
        list_to_write = []

        for i in range(1, 11):
            gene_rank_file = "%s_%s.txt" % (gene_rank_string, i)
            with open(gene_rank_file) as g_file:
                gene_rank_genes = g_file.read().splitlines()
            common_size = len(list(set(module_genes) & set(gene_rank_genes)))

            # take common_size amount of genes from this file and append them to a list
            n_random_genes = random.sample(gene_rank_genes, common_size)
            list_to_write.append(n_random_genes)

        flat_list = [item for sublist in list_to_write for item in sublist]
        # write the list to file as a representative randomized module
        with open('rand_%s_module_%s_prop.txt' % (out_string, j), 'w') as f:
            for one_gene in flat_list:
                f.write("%s\n" % one_gene)


if __name__ == "__main__":
    main()
