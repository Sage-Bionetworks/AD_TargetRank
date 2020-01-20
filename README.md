# AD_TargetRank
Rank and prioritize AD targets based on genetic and Genomic Evidence


## Repository Structure and Design:

All code is stored in ```code/```

The script to run the pipeline is ```code/Main.R``` it will draw-on and run the induvidual module scripts stored in ```code/modules/```

## Pipeline input:

Input to the pipeline is in YAML format. This format will specify the target list to load and run the pipeline on, as well as which modules scripts to run. It will also specify the synapse parentID and the folder name to create a synapse folder containing the results. Induvidual config files will be stored in ```configs/```. An example YAML input file will be displayed shortly.

## Output:

The compiled output will be a Wiki on the target synapse folder. This folder will also contain the output plots in induvidual files from the script. These will be stored in ```figures/``` and that entire folder is also specified in ```.gitignore``` to prevent this repo from growing too large.

## Setting up the pipeline run environment:

All dependencies will be contained within the docker image (idealy...) and can be run in an interactive RStudio session from an AWS instance utalizing the following commands:
```
sudo yum update -y
sudo amazon-linux-extras install docker
sudo service docker start
sudo usermod -aG docker <USER ID>
sudo docker run -v "/home/jgockley/AD_TargetRank:/home/jgockley/AD_TargetRank" -e USER=$(id -u ${USER}) -e PASSWORD=<PASSWORD> -d -p 8787:8787 <DockerImageID>

#Inside the docker image make the followning file
nano ~/.synapseConfig

#[authentication]
#username = <Synapse ID>
#password = <Synapse Password>

```

## Running the pipeline:

Sourcing the Inlitiazer with your yaml configuration file copies the path into the RMarkdown file and creates the Run.Rmd file. It then sources that file and knits the markdown file into the wiki of the synapse folder specified in the configuration script and copies the figures and tables into sub folders within the destination synapse folder specified in the configuration yaml file.

## Configuration file fields:

