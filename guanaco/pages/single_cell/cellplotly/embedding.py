import plotly.graph_objs as go
import pandas as pd
import numpy as np
import plotly.express as px


def plot_continuous_embedding(
    adata, embedding_key, color, x_axis=None, y_axis=None,
    transformation=None, order=None, color_map='Viridis',
    marker_size=5, opacity=1, annotation=None, axis_show=True
):
    """
    Plot a continuous feature (e.g., gene expression) on a 2D embedding.
    """
    embedding_prefixes = {
        'X_umap': 'UMAP', 'X_pca': 'PCA', 'X_tsne': 't-SNE',
        'X_diffmap': 'DiffMap', 'X_phate': 'PHATE', 'X_draw_graph_fa': 'FA'
    }
    embedding_prefix = embedding_prefixes.get(embedding_key, embedding_key.upper())
    embedding_data = adata.obsm[embedding_key]

    # Set column names for the embedding
    num_dimensions = embedding_data.shape[1]
    embedding_columns = [f'{embedding_prefix}{i + 1}' for i in range(num_dimensions)]
    embedding_df = pd.DataFrame(embedding_data, columns=embedding_columns)

    # Default x and y axis
    x_axis = x_axis or embedding_columns[0]
    y_axis = y_axis or (embedding_columns[1] if len(embedding_columns) > 1 else embedding_columns[0])

    # Extract gene expression and apply transformations
    gene_expression = adata[:, color].X.toarray().flatten()
    if transformation == 'log':
        gene_expression = np.log1p(gene_expression)
    elif transformation == 'z_score':
        gene_expression = (gene_expression - np.mean(gene_expression)) / np.std(gene_expression)

    embedding_df[color] = gene_expression

    # Sort order
    if order == 'max':
        embedding_df_sorted = embedding_df.sort_values(by=color)
    elif order == 'min':
        embedding_df_sorted = embedding_df.sort_values(by=color, ascending=False)
    elif order == 'random':
        embedding_df_sorted = embedding_df.sample(frac=1, random_state=315).reset_index(drop=True)
    else:
        embedding_df_sorted = embedding_df

    # Add annotation if provided
    if annotation:
        embedding_df_sorted[annotation] = adata.obs[annotation].values

    # Create scatter plot
    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=embedding_df_sorted[x_axis],
        y=embedding_df_sorted[y_axis],
        mode='markers',
        marker=dict(
            color=embedding_df_sorted[color],
            colorscale=color_map,
            cmin=embedding_df_sorted[color].min(),
            cmax=embedding_df_sorted[color].max(),
            size=marker_size,
            opacity=opacity,
            colorbar=dict(
                title=transformation if transformation else color,
                len = 0.8
            )
        ),
        text=embedding_df_sorted[annotation] if annotation else None,
        hovertemplate=f'{color}: %{{marker.color:.4f}}<br><extra></extra>'
    ))

    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        title=dict(
            text=f'<b>{color}</b>',
            x=0.5,
            y=0.9,
            xanchor='center'
        ),
        xaxis=dict(
            title=x_axis,
            showgrid=False,
            zeroline=False,
            scaleanchor='y',
            constrain='domain'
        ),
        yaxis=dict(
            title=y_axis,
            showgrid=False,
            zeroline=False,
            constrain='domain'
        )
    )

    fig.update_xaxes(
        showline=True, linewidth=2, linecolor='black',
        tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
    )
    fig.update_yaxes(
        showline=True, linewidth=2, linecolor='black',
        tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
    )

    return fig


# def plot_categorical_embedding(
#     adata, gene, embedding_key, color, x_axis=None, y_axis=None,
#     color_map=None, marker_size=5, opacity=1, legend_show='on legend',
#     axis_show=True
# ):
#     """
#     Plot a categorical variable (e.g., cell type) on a 2D embedding,
#     optionally showing a second variable (gene expression) on hover.
#     """
#     embedding_prefixes = {
#         'X_umap': 'UMAP', 'X_pca': 'PCA', 'X_tsne': 't-SNE',
#         'X_diffmap': 'DiffMap', 'X_phate': 'PHATE', 'X_draw_graph_fa': 'FA'
#     }
#     on_data = legend_show == 'on data'
#     embedding_data = adata.obsm[embedding_key]
#     prefix = embedding_prefixes.get(embedding_key, embedding_key.upper())
#     dims = [f'{prefix}{i + 1}' for i in range(embedding_data.shape[1])]

#     x_axis = x_axis or dims[0]
#     y_axis = y_axis or (dims[1] if len(dims) > 1 else dims[0])

#     df = pd.DataFrame(embedding_data, columns=dims)
#     df[color] = adata.obs[color].values
#     df[gene] = adata[:, gene].X.toarray().flatten()

#     unique_labels = sorted(df[color].unique())
#     color_map = color_map or px.colors.qualitative.Plotly
#     label_to_color = {label: color_map[i % len(color_map)] for i, label in enumerate(unique_labels)}

#     fig = go.Figure()

#     for label in unique_labels:
#         mask = df[color] == label
#         fig.add_trace(go.Scattergl(
#             x=df.loc[mask, x_axis],
#             y=df.loc[mask, y_axis],
#             mode='markers',
#             marker=dict(size=marker_size, color=label_to_color[label], opacity=opacity),
#             name=str(label),
#             customdata=np.stack([df.loc[mask, color], df.loc[mask, gene]], axis=-1),
#             hovertemplate=f'{color}: %{{customdata[0]}}<br>{gene}: %{{customdata[1]:.2f}}<extra></extra>',
#             showlegend=not on_data
#         ))

#     # Optional: add category labels on the plot
#     if on_data:
#         for label in unique_labels:
#             mask = df[color] == label
#             median_x = df.loc[mask, x_axis].median()
#             median_y = df.loc[mask, y_axis].median()
#             fig.add_annotation(
#                 x=median_x, y=median_y,
#                 text=f'<b>{label}</b>',
#                 showarrow=False,
#                 font=dict(size=12, color='black'),
#                 xanchor='center', yanchor='middle', opacity=0.9
#             )

#     fig.update_layout(
#         plot_bgcolor='white',
#         paper_bgcolor='white',
#         title=dict(
#             text=f'<b>{color}</b>',
#             x=0.5,
#             y=0.9,
#             xanchor='center'
#         ),
#         xaxis=dict(
#             title=x_axis,
#             showgrid=False,
#             zeroline=False,
#             scaleanchor='y',
#             constrain='domain',
#             tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
#         ),
#         yaxis=dict(
#             title=y_axis,
#             showgrid=False,
#             zeroline=False,
#             constrain='domain',
#             tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
#         )
#     )

#     fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
#     fig.update_yaxes(showline=True, linewidth=2, linecolor='black')

#     return fig
def plot_categorical_embedding(
    adata, gene, embedding_key, color,
    x_axis=None, y_axis=None,
    color_map=None, marker_size=5, opacity=1,
    legend_show='on legend', axis_show=True
):
    embedding_prefixes = {
        "X_umap": "UMAP", "X_pca": "PCA", "X_tsne": "t-SNE",
        "X_diffmap": "DiffMap", "X_phate": "PHATE", "X_draw_graph_fa": "FA"
    }
    on_data = legend_show == 'on data'

    # Prepare embedding coordinates
    embedding_data = adata.obsm[embedding_key]
    prefix = embedding_prefixes.get(embedding_key, embedding_key.upper())
    dims = [f"{prefix}{i+1}" for i in range(embedding_data.shape[1])]
    x_axis = x_axis or dims[0]
    y_axis = y_axis or (dims[1] if len(dims) > 1 else dims[0])

    # Prepare DataFrame
    df = pd.DataFrame(embedding_data, columns=dims)
    df[color] = adata.obs[color].values
    df[gene] = adata[:, gene].X.toarray().flatten()

    # Label & color mapping
    unique_labels = sorted(df[color].unique())
    color_map = color_map or px.colors.qualitative.Plotly
    label_to_color = {
        label: color_map[i % len(color_map)]
        for i, label in enumerate(unique_labels)
    }

    fig = go.Figure()

    # Add one trace per category
    for label in unique_labels:
        mask = df[color] == label
        fig.add_trace(go.Scattergl(
            x=df.loc[mask, x_axis],
            y=df.loc[mask, y_axis],
            mode='markers',
            marker=dict(
                size=marker_size,
                color=label_to_color[label],
                opacity=opacity,
            ),
            name=str(label),
            customdata=np.stack([df.loc[mask, color], df.loc[mask, gene]], axis=-1),
            hovertemplate=(
                f"{color}: %{{customdata[0]}}<br>"
                f"{gene}: %{{customdata[1]:.2f}}<extra></extra>"
            ),
            showlegend=not on_data
        ))

    # Add labels at cluster medians if requested
    if on_data:
        for label in unique_labels:
            mask = df[color] == label
            median_x = df.loc[mask, x_axis].median()
            median_y = df.loc[mask, y_axis].median()
            fig.add_annotation(
                x=median_x, y=median_y,
                text=f"<b>{label}</b>",
                showarrow=False,
                font=dict(size=12, color='black'),
                xanchor='center', yanchor='middle',
                opacity=0.9,
            )

    # Layout settings
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        title=dict(text=f"<b>{color}</b>", x=0.5, y=0.9, xanchor='center'),
        xaxis=dict(
            title=x_axis,
            showgrid=False, zeroline=False,
            scaleanchor='y', constrain='domain',
            tickfont=dict(color="rgba(0,0,0,0)" if not axis_show else "black")
        ),
        yaxis=dict(
            title=y_axis,
            showgrid=False, zeroline=False,
            constrain='domain',
            tickfont=dict(color="rgba(0,0,0,0)" if not axis_show else "black")
        ),
        legend=dict(
            orientation='v',
            itemsizing='constant',
            x=1.02, y=0.5,
            bgcolor='rgba(0,0,0,0)',
            itemclick='toggle',
            itemdoubleclick='toggleothers',
            font=dict(size=10)
        ) if not on_data else None,
        margin=dict(t=60, r=10, l=10, b=40)
    )

    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black')

    return fig

