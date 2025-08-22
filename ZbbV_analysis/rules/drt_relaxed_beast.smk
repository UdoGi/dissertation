import pandas as pd
import random
import xml.etree.ElementTree as et


# Parse XML
tree = et.parse(config['beast_xml'])
root = tree.getroot()
def yieldit():
    for taxon in root.findall('taxa')[0].findall('taxon'):
        taxon_id = taxon.attrib["id"]
        fasta_id = taxon_id.split(' ')[0]
        yield [fasta_id, fasta_id.split('_')[-2]]
            
df = pd.DataFrame(yieldit(), 
             columns=['fasta_id', 'year']
            )

# Parameters
N = 60  # number of replicates

sample_ids = list(df["fasta_id"])
dates = list(df["year"])

rule all:
    input:      
        expand("drt_relaxed/relaxed_randomized_{i}/relaxed_randomized_{i}.log", i=range(N))
        
rule shuffle_dates:
    output:
        xml="drt_relaxed/relaxed_randomized_{i}/relaxed_randomized_{i}.xml"
    run:
        # Shuffle dates
        randomized_dates = dates.copy()
        random.shuffle(randomized_dates)

        id_to_date = dict(zip(sample_ids, randomized_dates))

        # Find all <taxon> elements and update their <date> children
        for taxon in root.findall('taxa')[0].findall('taxon'):
            taxon_id = taxon.attrib["id"]
            fasta_id = taxon_id.split(' ')[0]
            if fasta_id in id_to_date:
                date_element = taxon.find("date")
                if date_element is not None:
                    date_element.set("value",str(id_to_date[fasta_id]))

        # Write out updated XML
        tree.write(output.xml)



rule run_beast:
    input:
        xml="drt_relaxed/relaxed_randomized_{i}/relaxed_randomized_{i}.xml"
    output:
        log="drt_relaxed/relaxed_randomized_{i}/relaxed_randomized_{i}.log"
    conda:
        "beast"
    threads: 
        1
    params:
        my_dir = lambda wildcards: f"drt_relaxed/relaxed_randomized_{wildcards.i}",
        i= lambda wildcards: wildcards.i
    shell:
        "cd {params.my_dir};"
        "beast relaxed_randomized_{params.i}.xml > relaxed_randomized_{params.i}.log;"
        "mv ancient_p_L_alignment_GTR_relaxed_tips.afa.log ancient_p_L_alignment_GTR_relaxed_tips_random_{params.i}.afa.log"
        