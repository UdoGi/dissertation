from Bio import Phylo
from pathlib import Path

map_d = {
    "NC_027717.1": "Hunter Island (NC_027717)",
    "Consensus_ALL_D10118_L_threshold_0_quality_20": 'Groningen 2018',
    "Consensus_ALL_D10758_L_threshold_0_quality_20": 'Sandau 1999',
    "Consensus_ALL_D10116_L_threshold_0_quality_20": "Hoogezand 2018",
    "Consensus_ALL_D10376_L_threshold_0_quality_20": 'Gusterath_2010',
    "Consensus_ALL_D9598_L_threshold_0_quality_20": 'Bad Lauterberg 2010',
    "Consensus_ALL_D10643_L_threshold_0_quality_20": 'Rossla 2011',
    "Consensus_ALL_D9593_L_threshold_0_quality_20": "Bad Salzdetfurth 2011",
    "Consensus_ALL_D9591_L_threshold_0_quality_20": "Diekholzen 2012",
    "Consensus_ALL_D6879_L_threshold_0_quality_20": "Penzlin 1919",
    "Consensus_ALL_D9595_L_threshold_0_quality_20": "Bad Münder 2012",
    "KF186494.1": "Malsoor (KF186494)",
    "KF186497.1": "Malsoor (KF186497)",
    "NC_021242.1": "Lone Star (NC_021242)",
    "JX961619.1": "Bhanja (JX961619)",
    "NC_022630.1": "Razdan (NC_022630)",
    "NC_027140.1": "Bhanja (NC_027140)",
    "NC_078070.1": "Kismayo (NC_078070)",
    "MN509914.2": "SFTS (MN509914)",
    "KR017835.1": "SFTS (KR017835)",
    "NC_043450.1": "SFTS (NC_043450.1)",
    "NC_043611.1": "Guertu (NC_043611)",
    "PP580190.1": "Kinna (PP580190)",
    "MZ440342.1": "Heartland (MZ440342)",
    "NC_078064.1": "Heartland (NC_078064)",
    "NC_078067.1": "Hunter Island (NC_078067.1)",
    "NC_014397.1": "Rift Valley Fever (NC_014397.1)",  # outgroup rift valley fever
    "MN823639.1": "Zwiesel bat bandavirus"
}

if __name__ == '__main__':

    with open(snakemake.log.stdout, "w") as stdout_log_file, open(snakemake.log.stderr, "w") as stderr_log_file:
        fn = Path(snakemake.input.tree)
        tree = Phylo.read(fn, format='newick')
        for clade in tree.get_terminals():  # Iterate over leaves
            clade.name = map_d[clade.name]
        Phylo.write(tree, snakemake.output.tree, 'newick')
