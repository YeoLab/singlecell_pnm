import matplotlib as mpl
import numpy as np
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
import seaborn as sns
from sklearn.decomposition import PCA

from anchor import MODALITY_TO_COLOR



def kmer_clustermap(kmer_zscores, retain=('included', 'bimodal', 'excluded'), 
                    row_filter=lambda x: x.var(axis=1) > (x.var(axis=1).mean() + 4*x.var(axis=1).std()),
                    **kwargs):
    data = kmer_zscores
    data = retain_cols(data, retain)
    print(data.shape)
    if row_filter is not None:
        data = data.loc[row_filter(data)]
    print(data.shape)

    intron_colors = make_intron_colors(data.columns)

    g = sns.clustermap(data.fillna(0), col_colors=intron_colors.values, method='ward', **kwargs)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0);
    return g


def make_intron_colors(columns):
    phenotype_colors = columns.map(lambda x: study.phenotype_to_color[x.split('_')[0]])
    modality_colors = columns.map(lambda x: MODALITY_TO_COLOR[x.split('_')[1]])
    direction_colors = columns.map(lambda x: direction_to_color[x.split('_')[-1].split('nt')[0].rstrip('0123456789')])
    intron_colors = pd.DataFrame([modality_colors, direction_colors, phenotype_colors], columns=columns)
    return intron_colors


def retain_cols(data, retain=('included', 'bimodal', 'excluded')):
    try:
        data = data[[col for col in data if any([r in col for r in retain])]]
    except TypeError:
        pass
    return data


def kmer_pcaplot(zscores, title, filename, retain=('included', 'excluded', 'bimodal'), transpose=False):
#     zscores = kmer_zscores[background]

    data = zscores.fillna(0)
    data = retain_cols(data, retain=retain)
    
    if transpose:
        data = data.T
        title += '-transposed'
        
    
    print(data.shape)
    

    pca = PCA(n_components=2)
    reduced = pd.DataFrame(pca.fit_transform(data), 
                           index=data.index)
    reduced.columns = reduced.columns.map(lambda x: 'component_{}'.format(x))
    # reduced = reduced
    print('\t', reduced.shape)
    reduced.index = reduced.index.map(lambda x: x.replace('T', 'U'))

    trace0 = go.Scatter(x=reduced.iloc[:, 0], y=reduced.iloc[:, 1], mode='markers', name='Motifs',
                        marker=dict(size=10, opacity=0.5, color='black'), text=reduced.index)

    
    lines = []
    if not transpose:
        metadata = pd.DataFrame(list(data.columns.map(lambda x: x.split('_'))), 
                                index=data.columns, columns=['phenotype', 'modality', 'location'])
        components = pd.DataFrame(pca.components_, columns=data.columns)
        print('\t', metadata.shape)
        scaling_factor = reduced.apply(np.linalg.norm, axis=1).max()
        for phenotype, phenotype_df in components.groupby(metadata['phenotype'], axis=1):
        #     linestyle = '-'
            if phenotype == 'iPSC':
                linestyle = 'solid'
            elif phenotype == 'NPC':
                linestyle = 'dash'
            else:
                linestyle = 'dot'
            for modality, modality_df in phenotype_df.groupby(metadata['modality'], axis=1):
                palette = map(mpl.colors.rgb2hex, reversed(sns.light_palette(MODALITY_TO_COLOR[modality], n_colors=3)))

                for color, (component, column) in zip(palette, modality_df.iteritems()):
                    x = [0, column[0]*scaling_factor]
                    y = [0, column[1]*scaling_factor]
                    lines.append(go.Scatter(x=x, y=y, mode='lines', name=component, text=component,
                                            line=dict(color=color, width=10, dash=linestyle)))

    plotly_data = [trace0] + lines

    layout = go.Layout(
    #     autosize=False, width=500, height=500,
    title=title,
    hovermode='closest',
    xaxis=dict(
        title='PC 1 explains {:d}% of variance'.format(int(pca.explained_variance_ratio_[0] * 100)),
        ticklen=5,
        zeroline=False,
        gridwidth=0,
        ),
    yaxis=dict(
        title='PC 2 explains {:d}% of variance'.format(int(pca.explained_variance_ratio_[1] * 100)),
        ticklen=5,
        zeroline=True,
        gridwidth=0,
        ),
    )

    fig = go.Figure(data=plotly_data, layout=layout)
    return py.iplot(fig, filename=filename)