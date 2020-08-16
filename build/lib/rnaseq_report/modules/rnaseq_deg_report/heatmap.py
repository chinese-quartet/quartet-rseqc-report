import plotly.graph_objects as go
# import plotly.figure_factory as ff
import rnaseq_report.modules.plotly_patch as pp

import numpy as np
from scipy.spatial.distance import pdist, squareform


def heatmap(df):

    df_array_s = np.array(df)
    df_array_g = df_array_s.transpose()

    # Initialize figure by creating upper dendrogram
    fig = pp.create_dendrogram(df_array_g, orientation='bottom')
    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'

    # Create Side Dendrogram
    dendro_side = pp.create_dendrogram(
        df_array_s, orientation='right')
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'

    # Add Side Dendrogram Data to Figure
    fig2 = fig
    for data in dendro_side['data']:
        fig2.add_trace(data)

    # Create Heatmap
    dendro_leaves_s = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves_g = fig['layout']['xaxis']['ticktext']
    dendro_leaves_s = list(map(int, dendro_leaves_s))
    dendro_leaves_g = list(map(int, dendro_leaves_g))

    heat_data = df_array_s[dendro_leaves_s, :]
    heat_data = heat_data[:, dendro_leaves_g]

    heatmap = [
        go.Heatmap(
            x=dendro_leaves_s,
            y=dendro_leaves_g,
            z=heat_data,
            colorscale='Blues'
        )
    ]

    heatmap[0]['x'] = dendro_side['layout']['yaxis']['tickvals']
    heatmap[0]['y'] = fig['layout']['xaxis']['tickvals']

    # Add Heatmap Data to Figure
    for data in heatmap:
        fig2.add_trace(data)

    # Edit Layout
    fig2.update_layout({'width': 800, 'height': 800,
                        'showlegend': False, 'hovermode': 'closest',
                        })

    # Edit xaxis
    fig2.update_layout(yaxis={'domain': [0, .85],
                              'mirror': False,
                              'showgrid': False,
                              'showline': False,
                              'zeroline': False,
                              'ticks': ""})
    # Edit xaxis2
    fig2.update_layout(yaxis2={'domain': [.825, .975],
                               'mirror': False,
                               'showgrid': False,
                               'showline': False,
                               'zeroline': False,
                               'showticklabels': False,
                               'ticks': ""})

    # Edit yaxis
    fig2.update_layout(xaxis={'domain': [.15, 1],
                              'mirror': False,
                              'showgrid': False,
                              'showline': False,
                              'zeroline': False,
                              'showticklabels': False,
                              'ticks': ""
                              })

    # Edit yaxis2
    fig2.update_layout(xaxis2={'domain': [0, .15],
                               'mirror': False,
                               'showgrid': False,
                               'showline': False,
                               'zeroline': False,
                               'showticklabels': False,
                               'ticks': ""})

    return fig2
