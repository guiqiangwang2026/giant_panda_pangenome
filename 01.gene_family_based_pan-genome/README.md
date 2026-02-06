### 1. Prepare gene cluster file

```shell
# Run Orthofinder to obtain the results of gene coding clustering
bash Orthofinder.sh
```

### 2. Gene family classification

```shell
# data preparation: Orthogroups.GeneCount.tsv + Orthogroups_UnassignedGenes.tsv
# copy the Orthogroups results file to the current path.
cp ../OrthoFinder/Results_Nov29/Orthogroups/Orthogroups.GeneCount.tsv ./
cp ../OrthoFinder/Results_Nov29/Orthogroups/Orthogroups_UnassignedGenes.tsv ./

# Run R script to generate pie charts and bar charts.
Rscript 00.Different_type_family.R
```

### 3. Pan and core genome size simulation 

```shell
# Preparing PanPG input files
# Run the R script to get PanGP.input
Rscript 01.prepare_PanGP_input.R

# Drawing using PanGP software
# download：https://sourceforge.net/projects/pangp/files/latest/download
# Running PanGP will generate PanGenomeData.txt
Rscript 02.fit_pan_core_plot.R
```

### 4. Analysis of the presence and absence (PAV) of genes

```shell
# Running R scripts for statistics and plotting
Rscript 04.PAV_heat_map.R
```
