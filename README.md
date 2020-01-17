# AD_TargetRank
Rank and prioritize AD targets based on genetic and Genomic Evidence


## Repository Structure and Design:

All code is stored in ```code/```

The script to run the pipeline is ```code/Main.R``` it will draw-on and run the induvidual module scripts stored in ```code/modules/```

## Pipeline input:

Input to the pipeline is in YAML format. This format will specify the target list to load and run the pipeline on, as well as which modules scripts to run. It will also specify the synapse parentID and the folder name to create a synapse folder containing the results. An example YAML input file will be displayed shortly.

## Output:

The compiled output will be a Wiki on the target synapse folder. This folder will also contain the output plots in induvidual files from the script. These will be stored in ```figures/``` and that entire folder is also specified in ```.gitignore```
