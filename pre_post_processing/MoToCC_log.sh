# MoToCC pre-processing

# to pre-process and save expression data for brain and kidney
# note example commands for dimensionality reduction, silhouette score, and percent composition plots within R scripts
Rscript exp_processing_brain.R
Rscript exp_processing_kidney.R


# compress expression, SNN, KNN
python3 compress_inputs.py norm_scale_exp_gene.csv SNN_graph.csv KNN_graph.csv out_string


# get cell ids, gene ids
# brain
cd brain/
head -n1 norm_scale_exp.csv | sed 's/,/\n/g' > gene_ids.txt
cut -d, -f1 RNA_graph.csv | tail -n+2 > cell_ids.txt

# kidney
# gene_ids must be gene symbols. 
# gene symbols download: http://uswest.ensembl.org/biomart/martview/c77c140f3051e795ceb0b839ecf55fbe 
# gene stable ID, gene name for hg38
cd kidney/
python3 header_matching.py \
norm_scale_exp.csv \
mart_export.txt \
norm_scale_header.csv
cat <(head -n1 norm_scale_header.csv) <(tail -n+2 norm_scale_exp.csv) > norm_scale_exp_gene.csv
head -n1 norm_scale_exp_gene.csv | sed 's/,/\n/g' > gene_ids.txt
cut -d, -f1 RNA_graph.csv | tail -n+2 > cell_ids.txt


# MoToCC
python3 MoToCC.py \
-e exp_norm_scale.hdf5 \
-sim sim_sparse.npz \
-knn knn_sparse.npz \
-m module.txt \
-c cell_ids.txt \
-g gene_ids.txt \
-k 100 \
-o outString

# return percent composition and percent capture 
# suppose you have solution files for multiple k 250 to 5,000 in steps of 250
for i in `seq 250 250 5000`
do
  sol="M2_${i}_LPonLSCC.txt"
  python3 get_labels.py $sol cell_metadata.csv OutFile_Composition
done

### Example post-processing for generation of representative randomized modules, calculation of naive s_h
# create randomized modules of same size as true module, with proportional expression levels 
python3 proportion_gene_exp.py ASD_plus_ID_M1_module.txt cortex_top_exp_genes M1
python3 proportion_gene_exp.py ASD_w_ID_M2_module.txt cortex_top_exp_genes M2
python3 proportion_gene_exp.py Optimal_SCN1A cortex_top_exp_genes SCN1A
python3 proportion_gene_exp.py ASD_mod13_Group1_valid.txt cortex_top_exp_genes ASD_mod13_Group1
python3 proportion_gene_exp.py Parkinson/C3_module_valid.txt Parkinson/parkinson_top_exp_genes_symbol C3 

# return naive s_h for selected cells < k 
python3 heuristic_summation.py codex/cell_metadata.csv ASD_plus_ID_M1_module.txt exp_norm_scale.hdf5 M1_Upregulate
python3 heuristic_summation.py codex/cell_metadata.csv ASD_plus_ID_M1_module.txt exp_norm_scale.hdf5 M1
python3 heuristic_summation.py codex/cell_metadata.csv  ASD_w_ID_M2_module.txt exp_norm_scale.hdf5 M2
python3 heuristic_summation.py codex/cell_metadata.csv  Optimal_SCN1A exp_norm_scale.hdf5 SCN1A
python3 heuristic_summation.py codex/cell_metadata.csv  ASD_mod13_Group1_valid.txt exp_norm_scale.hdf5 ASD_mod13_Group1
python3 heuristic_summation.py get_labels_input.txt C3_module_valid.txt  exp_norm_scale.hdf5 C3
