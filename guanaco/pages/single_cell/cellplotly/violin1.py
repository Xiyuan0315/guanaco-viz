from dash.exceptions import PreventUpdate
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
from guanaco.data_loader import color_config

def filter_data(adata, annotation, selected_labels):
    if selected_labels:
        # Filter the data based on the selected annotation labels
        cell_indices = adata.obs[annotation].isin(selected_labels)
        adata_filtered = adata[cell_indices]
    else:
        adata_filtered = adata
    return adata_filtered

def plot_violin1(adata, genes,labels, groupby, transformation = None, show_box=False, show_points=False , groupby_label_color_map = None):

    num_genes = len(genes)
    if num_genes ==0:
        raise PreventUpdate
    genes = [g for g in genes if g in adata.var_names]
    if not genes:
        return go.Figure()

    filtered_adata = filter_data(adata, groupby, labels)
    unique_labels = sorted(adata.obs[groupby].unique()) 
    if groupby_label_color_map is None:
        groupby_label_color_map = {
            label: color_config[i % len(color_config)] for i, label in enumerate(unique_labels)
        }

    fig = make_subplots(rows=num_genes, cols=1, shared_xaxes=True, 
                        vertical_spacing=0.02)  # No subplot titles as we add y-axis labels

    # List to collect annotations for gene names
    annotations = []
    y_positions = [(i / (num_genes-1))  for i in range(num_genes)] if num_genes >1 else [0.5]
    y_positions = y_positions[::-1]  # Reverse the list
    for i, gene in enumerate(genes):
        if gene in adata.var_names:
            gene_expression_matrix = filtered_adata[:, gene].X.toarray().flatten()
        else:
            raise ValueError(f"Gene '{gene}' not found in adata.var_names.")

        # transformation
        if transformation == 'log':
            gene_expression_matrix = np.log1p(gene_expression_matrix)  # log1p for log(1 + x) transformation
        elif transformation == 'z_score':
            gene_expression_matrix = (gene_expression_matrix - gene_expression_matrix.mean(axis=0)) / gene_expression_matrix.std(axis=0)

        df = pd.DataFrame({
            'Expression': gene_expression_matrix,
            groupby: filtered_adata.obs[groupby],
            'Gene': gene
        })

        points_mode = 'all' if show_points else False

        grouped_data = df.groupby(groupby,observed=True)
        # for only selected genes
        for label in labels:
            group_data = grouped_data.get_group(label)
            fig.add_trace(
                go.Violin(
                    y=group_data['Expression'],
                    x=group_data[groupby],
                    width = 0.8,
                    box_visible=show_box,
                    points=points_mode,
                    meanline_visible=True,
                    showlegend=(i == 0),  # Show legend only once
                    name=label,
                    spanmode='hard',
                    bandwidth=0.2,
                    fillcolor=groupby_label_color_map[label],
                    line_color='DarkSlateGrey',
                    hoveron='violins',
                    hoverinfo='y',
                    jitter=0.05,
                    scalemode='width',
                    scalegroup=label,
                    marker=dict(size=2),  # Set smaller dot size
                ),
                row=i + 1, col=1
            )
            fig.update_yaxes(range=[0, df['Expression'].max()], row=i + 1, col=1)

        # add gene name to the left
        annotations.append(dict(
            x=-0.1,
            y=y_positions[i],
            xref="paper",
            yref="paper",
            text=gene,
            showarrow=False,
            font=dict(size=10),
            xanchor='right' ))

    # Calculate appropriate height based on number of genes
    fig_height = max(400, num_genes * 120)  # 120px per gene, minimum 400px
    
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(size=10),
        showlegend=False, 
        annotations=annotations,
        height=fig_height,
        margin=dict(l=80, r=20, t=20, b=40)
    )

    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black')
    
    return fig