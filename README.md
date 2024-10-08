# PathExt_BRCA

**Project Overview**<br>
Breast cancers exhibit substantial transcriptional heterogeneity, posing a significant challenge to the prediction of treatment response and prognostication of outcomes. Especially, translation of TNBC subtypes to the clinic remains a work in progress, in part because of a lack of clear transcriptional signatures distinguishing the subtypes. Our recent network-based approach, PathExt, demonstrates that global transcriptional changes in a disease context are likely mediated by a small number of key genes, and these mediators may better reflect functional or translationally relevant heterogeneity. We apply PathExt to 1059 BRCA tumors and 112 healthy control samples across 4 subtypes to identify frequent, key-mediator genes in each BRCA subtype. We further implemented PathExt on Single cell transcriptomes of BRCA subtype tumors which revealed a subtype-specific distribution of PathExt-identified genes in multiple cell types from the tumor microenvironment. Next, Application of PathExt to a TNBC chemotherapy response dataset identified TNBC subtype-specific key genes and biological processes associated with resistance. We described putative drugs that target top novel genes potentially mediating drug resistance.

**##################### Instructions for the Users ###################**

In this study, we implemented netwrok based approach to identify differentially regulated paths (also referred as TopNets). These TopNets are of two types; (i) Activated; and (ii) Repressed. Further, from these TopNets, we identify top central genes.

PathExt requires node weights as an input for each gene in the corresponding sample. For more details, refer PathExt paper at "https://academic.oup.com/bioinformatics/article/37/9/1254/5952670?login=true"

**#################### Node Weight Computation ######################**

For computing node weights, user need to provide the input data in the csv file format. The input data should consist of 3 columns. First column should be gene list. Second column should be the value of those genes in control sample and Third column should be the value of those genes in the case. Order of the column is very important for the code.

It's recommended that user should provide the gene values in the quantile normalize form.

Next, run code **"Activated_node_weight.r"** to compute node weight for the Activated Network and **"Repressed_node_weight.r"** to compute node weight for the Repressed Network. The codes are provided in the **"code"** folder

The provided code compute value for one sample. User can run the code in loop for multiple samples.

**Code Usage:**

**/usr/local/bin/Rscript Activated_node_weight.r**      ##### For computing Activated Node Weight<br>
**/usr/local/bin/Rscript Repressed_node_weight.r**      ##### For computing Repressed Node Weight<br>

**########################### Running PathExt for Generating the TopNets ##########**

1. Compute the percentile threshold and q-score at which you user want minimum nubmer of nodes in the topnet. Run the following commands inside the folder where all the python codes are present.

**a. mkdir test_data/results/temp**<br>
**b. python3 node_weight_matrix_colname_Pijs.py test_data/input_data Sample1 test_data/human_PPIN.txt 0.1 2 1000 test_data/results/Activated_response test_data/results/temp/Pij**<br>
**c. python3 fdr_rand_pijs_boxcox.py test_data/results/temp test_data/results/Pij_zscores.txt**<br>
**d. rm -rf test_data/results/temp**<br>
**e. python3 try_different_thresholds_node_weight_matrix.py test_data/input_data Sample1 test_data/human_PPIN.txt 2 test_data/results/Pij_zscores.txt test_data/results/thresh_TopNet_sizes.txt**

2. After running the following commands, select the best values from the output file **"thresh_TopNet_sizes.txt"**. For example, user selected 0.01 as percentile and 0.05 as q-score. Now run the following commands for generating the topnets with user selected pecentile and q-score.

**mkdir test_data/results/temp**<br>

**python node_weight_matrix_colname_Pijs.py test_data/input_data Sample1 test_data/human_PPIN.txt 0.01 2 1000 test_data/results/Activated_response test_data/results/temp/Pij**<br>

**python fdr_rand_pijs_boxcox.py test_data/results/temp test_data/results/Pij_zscores.txt**<br>

**rm -rf test_data/results/temp**<br>

**python benjamini_hochberg_boxcox.py test_data/results/Pij_zscores.txt 0.05 test_data/results/Pij_zscores_fdr.txt**<br>

**python extract_fdr_network.py test_data/results/Activated_response test_data/results/Pij_zscores_fdr.txt 0.05 test_data/results/Activated_Response_TopNet.txt**

Here,

a) "input_data" is a tab seprated microarray data file with all the samples to be studied;<br>
b) "human_PPIN.txt" is the unweighted network file<br>
c) "Sample1" is the name of perturbation sample to study<br>
d) "0.01" is the percentile threshold<br>
e) "2" is the path length threshold<br>
f) "0.05" is the q-score cutoff<br>
g) "1000" is the number of randomizations<br>
h) "results" is the output directory<br>
i) "thresh_TopNet_sizes.txt" is the output file with all the percentile and q-score threshold.<br>
j) "Activated_Response" is the file name for base response network (we'll put it in the output directory)<br>
k) "Activated_Response_TopNet.txt" is the file name for TopNet (we'll put it in the output directory)<br>

**########################## Computing centrality score and top central genes ###########**

After generating the topnet file, compute the centrality score of each gene by running the following command

**python calc_ripple_centrality.py test_data/results/Activated_Response_TopNet.txt test_data/results/Activated_epicenter**

Next step is sorting the output file on the basis of "ripple_centrality" score and selecting the top required genes.<br>
We have provided the example ouput file in the **"test_data/results/"** folder for the reader.

<br>

**########################## Computing DEGs ###########**

DEGs in the current study were computed in 2 fashion.

**First Approach**<br>
For Upregulated genes, compute LogFC of Case/Control and select the genes with LogFC >2 and for Downregulated genes, select the genes with LogFC >-2.

**Second Approach**<br>
Computing DEGs using DESeq2, run the following command

**/usr/local/bin/Rscript DESeq2.r**
The input file comprises of raw reads for case and control and metadata file. These files are provided here by the name DESeq2_input and DESeq2_metadata. The output will comprise the output generated by the DESeq2 tool. 
The code and files are present in the code folder.

User should install DESeq2 library before running the command.

<br>

**############################ Creating Venn Diagram #######################################**
Create a comma separated file to generate venn diagram. Sample input file is provided "venn_input.csv"<br>
Code is present in the **"code"** folder

**Code usage:**

**/usr/local/bin/Rscript venn_diag.r**


**############################ Creating Boxplot #######################################**
Create a comma separated file to generate boxplot. Sample input file is provided "boxplot_input.csv"<br>
Code is present in the **"code"** folder

**Code usage:**

**/usr/local/bin/Rscript boxplot.r**


**############################ Creating Barplot #######################################**
Create a comma separated file to generate barplot. Sample input file is provided "barplot_input.csv"<br>
Code is present in the **"code"** folder

**Code usage:**

**/usr/local/bin/Rscript barplot.r**

**############################ Creating Heatmap #######################################**

Create a comma separated file to generate heatmap. Sample input file is provided "heatmap_input.csv" and the generated output is a jpg image "Heatmap.jpg"<br>
Code is present in the **"code"** folder

**Code usage:**

**/usr/local/bin/Rscript heatmap.r**

**############################ Creating DotPlot #######################################**

Create a comma separated file to generate Dotplot. Sample input file is provided "dotplot_input.csv" and the generated output is a jpg image "Dotplot.jpg"<br>
Code is present in the **"code"** folder

**Code usage:**

**/usr/local/bin/Rscript dotplot.r**

**############################ Creating Enrichment Plot among GO terms #######################################**

Clusterprofiler package was used for creating enrichment plots among GO terms. Code used is provided in the folder **"code"**
Run the command as

**/usr/local/bin/Rscript BP_clusterprofiler.r**<br>
Note that above code can be used to identify molecular functions. User need to change **ont="BP"** to **ont="MF"** <br>

**#################### Predicting Non-responder using SVM based model ######################**

SVM model was trained on GSE41998 dataset to predict non-responder to a given neoadjuvant chemotherapy using 21 gene signature expression as a feature. User can use this model to predict whether a patient will respond to the treatment or not (positive label classify as non-responder). Input file is provided by the name "ml_test.csv". Model is provided by the name **TNBC_finalized_model.sav**. All the 3 files (code, test file and model is present in the maiin directory)

**Code usage:**
**python SVM_predict.py**


**########################## Virtual Screening #########################################**

For virtual screening, we downloaded the 3D structure of the protein from RCSB-PDB database and ligands in SMILES (.smi) format from ZINC database.<br>
We used AutoDock Vina for virtual screening experiments. The software requires the ligand file in .mol2 file format. Following command was used to convert the ligand from .smi to .mol2 file format using obabel software.

**obabel -i smi mol1 -o mol2 --gen3D -O test.mol2**

Once the ligand was converted into .mol2 file format, we prepared the protein for docking by the following command. This command will select the required chain for docking, and will remove the heteroatoms and unwanted metals from the structure.

**pdb_selchain -A 4qcl.pdb | pdb_delhetatm | pdb_tidy |grep "^ATOM" > 4qcl_processed.pdb**

Once the protein and ligands are ready, we need to convert them for docking. User should download the Autodock Vina software and install it.
Now we convert the receptor and ligand files into the .pdbqt file format. This file format is used by Vina for docking. Command for the same is

**python /usr/local/apps/Autodock/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -r receptor.pdb -o receptor.pdbqt**<br>
**python /usr/local/apps/Autodock/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l ligand.mol2 -o ligand.pdbqt**

After converting it into .pdbwt file format, run the code vina.pl as given below

**vina --config conf.txt --ligand ligand.pdbqt --log ligand_log.out**

This code requires the **"conf.txt file"** where we define the grid size of the box and the coordinates of the protein where we want to dock the ligand. This can be done using UCSF chimera tool.

**Tools and packages used:**<br>
Python 3.6.9<br>
R 3.6<br>
Pandas 0.25.3<br>
Networkx 1.11<br>
Numpy 1.17.4<br>
Statsmodel<br>
Random<br>
Sys<br>
Math<br>
ggplot2<br>
dplyr<br>
tidyverse<br>
pheatmap<br>
loess<br>
Clusterprofiler<br>
openbabel<br>
AutoDock Vina<br>


