
## Code related to the analysis of the ZbbV genomes

### Prerequisites
1. Install the necessary requirements in _../envs/all.yaml_ into a new conda environment  
`conda env create --name {your_chosen_environment_name} --file envs/all.yaml`
2. Download [IQ-TREE](http://www.iqtree.org/) and adjust the path to the binary in the config.yaml


### Running the snakemake pipeline
1. adapt the **config.yaml** file so that it points to the desired locations on your filesystem
2. navigate to the _rules_ director and execute the snakefile **full_pipeline.smk**  
`snakemake -s full_pipeline.smk --cores 8 --configfile ../config.yaml --use-conda`

## Snakemake rules overview
### merge_fastq_per_case.smk
Merges all sequencing data of an individual animal. This can be different sequencing runs of the same organ or different organs of the same animal

### preparations.smk
creates and downloads reference sequences that are subsequently used to remove unwanted reads from the fastqs.

### create_ancient_reference.smk
creates the genome reference of the ancient sample via de-novo assembly, which was the basis for all subsequent work

### create_recent_consensus.smk
Maps the sequence data of the recent samples against the ancient sample and calls a consensus

### add_sanger_data.smk
For some samples, additional sanger sequencing was done to fill gaps. This data is incorporated here.

### rename_consensi_MSA.smk
Renames the CDS of each sample using the samples location and collection year.
For each segment an alignment using muscle is then performed and iqtree is called on 
each of the alignments to build an ML tree.

### create_overview_tree.smk
Downloads references of bandaviruses, aligns with muscle and builds a tree with iqtree.
Renames the tip labels and runs IQtree again to compute alrt-bootstrap values using the best
selected model from the run before.

### full_pipeline.smk
This Snakefile includes all the rules above and executes them consecutively.
