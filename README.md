# Instructions on how to parse the Human CellAtlas

```
Rscript download.R /vast/scratch/users/mangiola.s/human_cell_atlas

```
# Execute makeflow

```


conda activate cctools-env
makeflow -J 200 -T slurm get_metadata.makeflow 

``