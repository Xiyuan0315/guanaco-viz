import numpy as np
import plotly.graph_objs as go
from dash.exceptions import PreventUpdate
import pandas as pd

def plot_dot_matrix(adata, genes, groupby, selected_labels, aggregation='mean', transformation=None,standardization = None, vmin=None, vmax=None, expression_threshold=0, color_map='Viridis', plot_type='dotplot'):
    valid_genes = [gene for gene in genes if gene in adata.var_names]
    if not valid_genes:
        raise PreventUpdate

    # Filter data based on selected labels
    if selected_labels:
        cell_indices = adata.obs[groupby].isin(selected_labels)
        adata = adata[cell_indices]

    expression_data = adata.to_df()[valid_genes].copy()
    expression_data[groupby] = adata.obs[groupby].values

    if transformation == 'log':
        expression_data[valid_genes] = np.log1p(expression_data[valid_genes])
    elif transformation == 'z_score':
        expression_data[valid_genes] = (expression_data[valid_genes] - expression_data[valid_genes].mean()) / expression_data[valid_genes].std()

    aggregated_data = expression_data.groupby(groupby, observed=True).agg(aggregation)
        
    if standardization == 'var':
        # Standardize per gene (each column)
        aggregated_data = (aggregated_data - aggregated_data.min()) / (aggregated_data.max() - aggregated_data.min())
    elif standardization == 'group':
        # Standardize per group (each row)
        aggregated_data = (aggregated_data.sub(aggregated_data.mean(axis=1), axis=0)
                                        .div(aggregated_data.std(axis=1), axis=0))
    else:
        pass

    expression_binary = expression_data[valid_genes] > expression_threshold
    fraction_expressing = expression_binary.groupby(expression_data[groupby], observed=True).mean()

    # Use actual groups present in the data after filtering
    if selected_labels:
        # Keep the order of selected_labels but only include those actually present
        groups = [label for label in reversed(selected_labels) if label in aggregated_data.index]
    else:
        groups = list(reversed(aggregated_data.index))
    
    vmin = vmin if vmin is not None else aggregated_data[valid_genes].min().min()
    vmax = vmax if vmax is not None else aggregated_data[valid_genes].max().max()

    if plot_type == 'dotplot':
        df_expression = aggregated_data.reset_index().melt(id_vars=[groupby], value_vars=valid_genes, var_name='gene', value_name='expression')
        df_fraction = fraction_expressing.reset_index().melt(id_vars=[groupby], value_vars=valid_genes, var_name='gene', value_name='fraction')
        df_merged = pd.merge(df_expression, df_fraction, on=[groupby, 'gene'])

        max_fraction = df_merged['fraction'].max()
        scale_params = (
            (150, [0.1, 0.075, 0.05, 0.025]) if max_fraction < 0.1 else
            (90, [0.2, 0.15, 0.1, 0.05]) if max_fraction < 0.2 else
            (40, [0.4, 0.3, 0.2, 0.1]) if max_fraction < 0.4 else
            (35, [0.5, 0.4, 0.3, 0.2]) if max_fraction < 0.5 else
            (30, [0.6, 0.5, 0.4, 0.3]) if max_fraction < 0.6 else
            (25, [0.7, 0.6, 0.5, 0.4]) if max_fraction < 0.7 else
            (20, [0.8, 0.6, 0.3, 0.1]) if max_fraction < 0.8 else
            (17, [1.0, 0.75, 0.5, 0.25])
        )
        scale, size_legend_values = scale_params

        marker_sizes = df_merged['fraction'] * scale
        size_legend_sizes = [s * scale for s in size_legend_values]
        custom_data = np.stack([df_merged['expression'], df_merged['fraction']], axis=-1)

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=df_merged['gene'],
            y=df_merged[groupby],
            mode='markers',
            showlegend=False,
            marker=dict(
                size=marker_sizes,
                color=df_merged['expression'],
                colorscale=color_map,
                cmin=vmin,
                cmax=vmax,
                colorbar=dict(
                    title=f'{aggregation.capitalize()} Expression ({transformation})' if transformation != 'None' else f'{aggregation.capitalize()} Expression',
                    tickfont=dict(color='DarkSlateGrey', size=10),
                    len=0.5,
                    yanchor="bottom",
                    y=0
                )
            ),
            line=dict(color='black', width=0.5),
            customdata=custom_data,
            hovertemplate=(
                'Gene: %{x}<br>'
                f'{groupby}: %{{y}}<br>'
                'Expression: %{customdata[0]:.4f}<br>'
                'Fraction: %{customdata[1]:.4f}<extra></extra>'
            )
        ))

        fig.add_trace(go.Scatter(
            x=[None], y=[""],
            mode="text",
            text=["Fraction of Cells"],
            textposition="top center",
            textfont=dict(size=10, color="black"),
            showlegend=False
        ))

        for size, value in zip(size_legend_sizes, size_legend_values):
            fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(size=size, color="grey", opacity=1),
                name=f"{'Percent Expressed <br>' if size == size_legend_sizes[0] else ''}{value:.2f}",
                legendgroup="size_legend",
                showlegend=True
            ))

        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            xaxis=dict(showline=True, linewidth=2, linecolor='black', showgrid=False),
            yaxis=dict(showline=True, linewidth=2, linecolor='black', showgrid=False, categoryorder='array', categoryarray=groups)
        )

    else:
        fig = go.Figure(data=go.Heatmap(
            z=aggregated_data[valid_genes].values,
            x=valid_genes,
            y=aggregated_data.index,
            colorscale=color_map,
            zmin=vmin,
            zmax=vmax,
            colorbar=dict(
                title=dict(text=f'{aggregation.capitalize()} Expression<br>({transformation})', font=dict(size=14)),
                tickfont=dict(size=14),
                outlinewidth=1
            ),
            showscale=True,
            xgap=0.8,
            ygap=0.8,
            hovertemplate=(
                f'{groupby}: %{{y}}<br>'
                'Gene: %{x}<br>'
                'Expression Level: %{z:.4f}<br><extra></extra>'
            )
        ))

        fig.update_layout(
            title=dict(x=0.5, xanchor='center', font=dict(size=16, color='DarkSlateGrey')),
            plot_bgcolor='white',
            paper_bgcolor='white',
            font=dict(size=14),
            showlegend=False,
            yaxis=dict(categoryorder='array', categoryarray=groups)
        )

    return fig
