## Code related to the discovery pipeline used in my dissertation

### Prerequisites
1. Install the necessary requirements in _../envs/all.yaml_ into a new conda environment  
`conda env create --name {your_chosen_environment_name} --file envs/all.yaml`


### Running the snakemake pipeline
1. look up the path to directory the FASTQ files are located
2. execute the snakefile **full_pipeline.smk** and specify the path to the FASTQ files  
`snakemake --cores 8 --config path_to_samples={FASTQ_directory} --use-conda`
