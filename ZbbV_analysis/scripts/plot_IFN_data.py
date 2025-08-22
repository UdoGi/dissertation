import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import sys


def plot_activation(df, index, colors, fn):
    sns.set(font_scale=2.6)
    fig, ax = plt.subplots(figsize=(20, 10))
    sns.barplot(data=df.loc[index].reset_index(), x="my_index", y="% activation", ax=ax)
    for bar, color in zip(ax.patches, colors):
        bar.set_color(color)
    _ = plt.xticks(rotation=90, fontsize=12)
    plt.savefig(fn, format="svg")


def plot_figure4(df, fn_out):
    cols_fig4 = [
        ("-", 'empty vector', 100),
        ("+", 'empty vector', 300),
        ('+', '3xFLAG_SFTSV NSs', 30),
        ('+', '3xFLAG_SFTSV NSs', 100),
        ('+', '3xFLAG_SFTSV NSs', 300),
        ('+', 'Ancient-NSs', 30),
        ('+', 'Ancient-NSs', 100),
        ('+', 'Ancient-NSs', 300),
        ('+', 'Recent-NSs', 30),
        ('+', 'Recent-NSs', 100),
        ('+', 'Recent-NSs', 300),
        ('+', '3xFLAG_HALO', 100),
    ]
    colors_fig4 = ['grey'] * 2 + ['orange'] * 3 + \
                  [(0.617762399077278, 0.6021376393694733, 0.7834525182622069)] * 3 + \
                  [(0.3568166089965398, 0.20525951557093425, 0.5856978085351787)] * 3 + \
                  [(0.9882352941176471, 0.6261437908496732, 0.5084967320261438)]

    plot_activation(df, cols_fig4, colors_fig4, fn_out)


def plot_supp_figure(df, fn_out, is_recent=False):
    cols_ancient = [
        ("-", '3xFLAG_HALO', 100),
        ("-", '3xFLAG_SFTSV NSs', 100),
        ("-", 'Ancient-NSs', 100),
        ("-", '3xFLAG Ancient-NSs', 100),
        ("-", 'Ancient-NSs-3xFLAG', 100),
        ("-", 'empty vector', 100),
        ("+", 'empty vector', 100), ("+", 'empty vector', 300),
        ("+", '3xFLAG_HALO', 30), ("+", '3xFLAG_HALO', 100), ("+", '3xFLAG_HALO', 300),
        ("+", '3xFLAG_SFTSV NSs', 30), ("+", '3xFLAG_SFTSV NSs', 100), ("+", '3xFLAG_SFTSV NSs', 300),
        ("+", 'Ancient-NSs', 30), ("+", 'Ancient-NSs', 100), ("+", 'Ancient-NSs', 300),
        ("+", '3xFLAG Ancient-NSs', 30), ("+", '3xFLAG Ancient-NSs', 100), ("+", '3xFLAG Ancient-NSs', 300),
        ("+", 'Ancient-NSs-3xFLAG', 30), ("+", 'Ancient-NSs-3xFLAG', 100), ("+", 'Ancient-NSs-3xFLAG', 300),
    ]
    colors_full = ['grey'] * 6 + [(0.9882352941176471, 0.6261437908496732, 0.5084967320261438)] * 2 + \
                  [(0.21568627450980393, 0.5294117647058824, 0.7542483660130719)] * 3 + ['orange'] * 3 + \
                  [(0.617762399077278, 0.6021376393694733, 0.7834525182622069)] * 3 + \
                  [(0.47320261437908495, 0.43267973856209146, 0.6993464052287581)] * 3 + \
                  [(0.3568166089965398, 0.20525951557093425, 0.5856978085351787)] * 3

    def replace_ancient(plasmid):
        return plasmid.replace('Ancient', 'Recent') if 'Ancient' in plasmid else plasmid

    if is_recent:
        cols_recent = [(_, replace_ancient(plasmid), conc) for _, plasmid, conc in cols_ancient]
        plot_activation(df, cols_recent, colors_full, fn_out)
    else:
        plot_activation(df, cols_ancient, colors_full, fn_out)


def prepare_dataframe(fn_in):
    df = pd.read_csv(fn_in)
    # dividing firefly luciferase bioluminescence value by background renilla luciferase bioluminescence value
    df['raw_ffluc_rluc_ratio'] = df.ffluc / df.rluc
    # create a new name for an experiment that ignores the suffix indicating if the experiment on that day was
    # with ancient or recent NSs. This allows to merge the redundant control measurements.
    df['experiment_together'] = df.experiment.apply(lambda x: x.split('_')[0])

    # aggregating on the level of experiment, s.t.
    # the two (ancient/recent) controls and sftsv are merged
    dfp = df.groupby(["MAVS", 'expression_plasmid',
                      "Amount of plamid (ng)", 'experiment_together'])[
        ['ffluc', 'rluc', 'raw_ffluc_rluc_ratio']
    ].mean().reset_index() \
        .set_index(["MAVS", 'expression_plasmid',
                    "Amount of plamid (ng)"])

    # all measurements are normalized to the activation level of the empty vector when MAVS is overexpressed
    # the created dictionary gives for each experiment this normalization value of the empty vector.
    empty_vec_pos_lookup = dict(dfp.loc[("+", 'empty vector', 300)][
                                    ['experiment_together', 'raw_ffluc_rluc_ratio']].values)

    def compute_perc_act(row):
        return 100 * row.raw_ffluc_rluc_ratio / empty_vec_pos_lookup[row.experiment_together]

    dfp['% activation'] = dfp.apply(compute_perc_act, axis=1)
    dfp['my_index'] = ["\n".join(map(str, x)).replace('NSs', '').replace('-', '')
                       for x in dfp.index]

    return dfp


if __name__ == "__main__":
    sys.stderr = open(snakemake.log.stderr, 'w')
    sys.stdout = open(snakemake.log.stdout, 'w')

    df = prepare_dataframe(snakemake.input.csv)
    plot_figure4(df, fn_out=snakemake.output.fn_fig4)
    plot_supp_figure(df, snakemake.output.fn_suppFig1, is_recent=False)
    plot_supp_figure(df, snakemake.output.fn_suppFig2, is_recent=True)

