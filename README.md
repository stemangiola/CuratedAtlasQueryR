# Instructions on how to parse the Human CellAtlas

# Step 1 download

```
Rscript download.R /vast/scratch/users/mangiola.s/human_cell_atlas

```
# Get gene names

```
Rscript get_gene_names.R /vast/scratch/users/mangiola.s/human_cell_atlas

```

# Create metadata

```

conda activate cctools-env
makeflow -J 200 -T slurm get_metadata.makeflow 

```

# Step 3 split files

```

conda activate cctools-env
makeflow -J 200 -T slurm split_files.makeflow 

```

# Step 4 consolidate files

```

conda activate cctools-env
makeflow -J 200 -T slurm light_files.makeflow 

```