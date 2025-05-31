from dash import dcc, html, Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import json
from pathlib import Path

from guanaco.pages.single_cell.cellplotly.heatmap1 import plot_heatmap1
from guanaco.pages.single_cell.cellplotly.heatmap2 import plot_heatmap2
from guanaco.pages.single_cell.cellplotly.violin1 import plot_violin1
from guanaco.pages.single_cell.cellplotly.violin2 import plot_violin2
from guanaco.pages.single_cell.cellplotly.stacked_bar import plot_stacked_bar
from guanaco.pages.single_cell.cellplotly.dotmatrix import plot_dot_matrix

# Load color palettes
# Load cvd_color.json from the same directory as this module
cvd_color_path = Path(__file__).parent / "cvd_color.json"
with open(cvd_color_path, "r") as f:
    palette_json = json.load(f)
palette_names = list(palette_json["color_palettes"].keys())

from guanaco.pages.single_cell.mod021_heatmap import generate_heatmap_layout
from guanaco.pages.single_cell.mod022_violin import generate_violin_layout
from guanaco.pages.single_cell.mod023_dotplot import generate_dotplot_layout
from guanaco.pages.single_cell.mod024_stacked_bar import generate_stacked_bar_layout

import warnings
warnings.filterwarnings('ignore', message='.*observed=False.*')


def filter_data(adata, annotation, selected_labels):
    if selected_labels:
        cell_indices = adata.obs[annotation].isin(selected_labels)
        adata_filtered = adata[cell_indices]
    else:
        adata_filtered = adata
    return adata_filtered


def generate_left_control(default_gene_markers, label_list, prefix):
    genes_selection = dcc.Dropdown(
        id=f'{prefix}-single-cell-genes-selection',
        options=[{'label': gene, 'value': gene} for gene in default_gene_markers],
        value=default_gene_markers,
        multi=True,
        style={'marginBottom': '15px', 'font-size': '12px'},
        className='custom-dropdown'
    )

    annotation_filter = dcc.Dropdown(
        id=f'{prefix}-single-cell-annotation-dropdown',
        options=[{'label': label, 'value': label} for label in label_list],
        value=label_list[len(label_list) // 2],
        style={'marginBottom': '15px'},
        clearable=False,
        className='custom-dropdown'
    )

    label_list_selection = dcc.Dropdown(
        id=f'{prefix}-single-cell-label-selection',
        multi=True,
        style={'marginBottom': '15px', 'font-size': '12px'},
        className='custom-dropdown'
    )

    # Discrete color map dropdown
    discrete_color_map_dropdown = dcc.Dropdown(
        id=f"{prefix}-discrete-color-map-dropdown",
        options=[{"label": name, "value": name} for name in palette_names],
        value=None,
        placeholder="Default color",
        style={"marginBottom": "15px"},
        clearable=True,
        className="custom-dropdown"
    )

    return html.Div([
        html.Label('Select Variables:', style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        genes_selection,
        html.Label('Select Annotation:', style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        annotation_filter,
        html.Label('Select Labels:', style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        label_list_selection,
        html.Label('Discrete ColorMap:', style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        discrete_color_map_dropdown
    ])


def generate_single_cell_tabs(adata, default_gene_markers, discrete_label_list, prefix):

    tabs = dcc.Tabs([
        dcc.Tab(label='Heatmap', value='heatmap-tab', children=[
            html.Div(generate_heatmap_layout(adata, prefix))
        ]),
        dcc.Tab(label='Violin Plot', value='violin-tab', children=[
            html.Div(generate_violin_layout(adata, default_gene_markers, discrete_label_list, prefix))
        ]),
        dcc.Tab(label='Dotplot', value='dotplot-tab', children=[
            html.Div(generate_dotplot_layout(prefix))
        ]),
        dcc.Tab(label='Stacked Bar', value='stacked-bar-tab', children=[
            html.Div(generate_stacked_bar_layout(discrete_label_list,prefix))
        ]),
    ], id=f'{prefix}-single-cell-tabs', value='heatmap-tab', className='custom-tabs')

    return dbc.Row([
        dbc.Col(
            generate_left_control(default_gene_markers, discrete_label_list, prefix),
            xs=12, sm=12, md=4, lg=4, xl=2
        ),
        dbc.Col(tabs, xs=12, sm=12, md=8, lg=8, xl=10)
    ], style={'marginBottom': '50px'})


def maker_vis_callback(app, adata, prefix):
    @app.callback(
        Output(f'{prefix}-violin2-group-selection', 'options'),
        Output(f'{prefix}-violin2-group-selection', 'value'),
        Input(f'{prefix}-multi-class-selection', 'value')
    )
    def update_group_labels(selected_column):
        unique_labels = sorted(adata.obs[selected_column].dropna().unique(), key=str)
        options = [{'label': str(label), 'value': str(label)} for label in unique_labels]
        values = [str(label) for label in unique_labels]
        return options, values

    @app.callback(
        Output(f'{prefix}-single-cell-genes-selection', 'options'),
        Input(f'{prefix}-single-cell-genes-selection', 'search_value'),
        State(f'{prefix}-single-cell-genes-selection', 'value')
    )
    def update_genes_dropdown(search_value, value):
        if not search_value:
            raise PreventUpdate
        label_list = adata.var_names.to_list()
        matching_labels = [label for label in label_list if search_value.lower() in label.lower()]
        selected_labels = value if value else []
        all_labels = list(set(selected_labels + matching_labels[:10]))
        return [{'label': label, 'value': label} for label in all_labels]

    @app.callback(
        [Output(f'{prefix}-single-cell-label-selection', 'options'),
         Output(f'{prefix}-single-cell-label-selection', 'value')],
        Input(f'{prefix}-single-cell-annotation-dropdown', 'value')
    )
    def update_labels_based_on_annotation(selected_annotation):
        unique_labels = adata.obs[selected_annotation].unique().tolist()
        label_options = [{'label': label, 'value': label} for label in unique_labels]
        label_values = sorted(unique_labels)
        return label_options, label_values

    @app.callback(
        Output(f'{prefix}-heatmap', 'figure'),
        [Input(f'{prefix}-single-cell-genes-selection', 'value'),
         Input(f'{prefix}-single-cell-annotation-dropdown', 'value'),
         Input(f'{prefix}-single-cell-label-selection', 'value'),
         Input(f'{prefix}-heatmap-transformation', 'value'),
         Input(f'{prefix}-heatmap-colorscale-dropdown', 'value'),
         Input(f'{prefix}-heatmap-label-dropdown', 'value'),
         Input(f'{prefix}-discrete-color-map-dropdown', 'value')]
    )
    def update_heatmap(selected_genes, selected_annotation, selected_labels, transformation, heatmap_color, secondary_annotation, discrete_color_map):
        filtered_adata = filter_data(adata, selected_annotation, selected_labels)
        
        # Use heatmap2 if secondary annotation is different from primary annotation and not "None"
        if secondary_annotation and secondary_annotation != 'None' and secondary_annotation != selected_annotation:
            # Create color maps for consistent coloring
            unique_labels1 = sorted(adata.obs[selected_annotation].unique())
            unique_labels2 = sorted(adata.obs[secondary_annotation].unique())
            
            # Use discrete color map if selected, otherwise use default
            if discrete_color_map:
                discrete_palette = palette_json["color_palettes"][discrete_color_map]
                groupby1_label_color_map = {
                    label: discrete_palette[i % len(discrete_palette)] for i, label in enumerate(unique_labels1)
                }
                groupby2_label_color_map = {
                    label: discrete_palette[i % len(discrete_palette)] for i, label in enumerate(unique_labels2)
                }
            else:
                from guanaco.data_loader import color_config
                groupby1_label_color_map = {
                    label: color_config[i % len(color_config)] for i, label in enumerate(unique_labels1)
                }
                groupby2_label_color_map = {
                    label: color_config[i % len(color_config)] for i, label in enumerate(unique_labels2)
                }
            
            return plot_heatmap2(
                adata=filtered_adata,
                genes=selected_genes,
                groupby1=selected_annotation,
                groupby2=secondary_annotation,
                labels=selected_labels,
                log=(transformation == 'log'),
                z_score=(transformation == 'zscore'),
                boundary=1,
                color_map=heatmap_color,
                groupby1_label_color_map=groupby1_label_color_map,
                groupby2_label_color_map=groupby2_label_color_map
            )
        else:
            # Use heatmap1 for single annotation
            # Create color map for annotation bar if discrete color map is selected
            groupby_label_color_map = None
            if discrete_color_map:
                discrete_palette = palette_json["color_palettes"][discrete_color_map]
                unique_labels = sorted(adata.obs[selected_annotation].unique())
                groupby_label_color_map = {
                    label: discrete_palette[i % len(discrete_palette)] for i, label in enumerate(unique_labels)
                }
            
            return plot_heatmap1(
                adata=filtered_adata,
                genes=selected_genes,
                labels=selected_labels,
                adata_obs=adata.obs,
                groupby=selected_annotation,
                transformation=transformation,
                boundary=1,
                color_map=heatmap_color,
                groupby_label_color_map=groupby_label_color_map
            )

    @app.callback(
        Output(f'{prefix}-violin-plot1', 'figure'),
        [Input(f'{prefix}-single-cell-tabs', 'value'),
         Input(f'{prefix}-single-cell-genes-selection', 'value'),
         Input(f'{prefix}-single-cell-annotation-dropdown', 'value'),
         Input(f'{prefix}-single-cell-label-selection', 'value'),
         Input(f'{prefix}-violin-log-or-zscore', 'value'),
        Input(f'{prefix}-show-box1', 'value'),
         Input(f'{prefix}-show-scatter1', 'value'),
         Input(f'{prefix}-discrete-color-map-dropdown', 'value')]
    )
    def update_violin1(selected_tab, selected_genes, selected_annotation, 
                       selected_labels, transformation, show_box_plot, show_scatter1, discrete_color_map):
        if selected_tab == 'violin-tab':
            # Use discrete color map if selected, otherwise use default
            color_map = None
            if discrete_color_map:
                discrete_palette = palette_json["color_palettes"][discrete_color_map]
                unique_labels = sorted(adata.obs[selected_annotation].unique())
                color_map = {
                    label: discrete_palette[i % len(discrete_palette)] for i, label in enumerate(unique_labels)
                }
            
            fig = plot_violin1(
                adata,
                selected_genes,
                selected_labels,
                groupby=selected_annotation,
                transformation=transformation,
                show_box='show' in show_box_plot if show_box_plot else False,
                show_points='show' in show_scatter1 if show_scatter1 else False,
                groupby_label_color_map=color_map
            )
            num_genes = len(selected_genes)
            num_categories = len(selected_labels)
            fig.update_layout(
                height=min(1000, max(300, 80 * num_genes)),
                width=min(500, max(200, 110 * num_categories)),
                margin=dict(l=130, r=10, t=30, b=30)
            )
            return fig
        return go.Figure()

    @app.callback(
        Output(f'{prefix}-p-value-method', 'options'),
        Output(f'{prefix}-p-value-method', 'value'),
        Input(f'{prefix}-p-value-selection', 'value')
    )
    def update_p_value_method(selection):
        if selection == 'multi':
            options = [
                {'label': 'Kruskal-Wallis', 'value': 'kw-test'},
                {'label': 'ANOVA', 'value': 'anova'},
                {'label': 'Two-way ANOVA', 'value': 'two-way-anova'}
            ]
            value = 'kw-test'
        elif selection == 'binary':
            options = [
                {'label': 'Mann-Whitney U Test', 'value': 'mwu-test'},
                {'label': 'T-test', 'value': 'ttest'}
            ]
            value = 'mwu-test'
        else:
            options = [{'label': 'None', 'value': 'none'}]
            value = 'none'
        return options, value

    @app.callback(
        Output(f'{prefix}-violin2-gene-selection', 'options'),
        Input(f'{prefix}-violin2-gene-selection', 'search_value')
    )
    def update_violin_genes_dropdown(search_value):
        if not search_value:
            raise PreventUpdate
        label_list = adata.var_names.to_list()
        matching_labels = [label for label in label_list if search_value.lower() in label.lower()]
        return [{'label': label, 'value': label} for label in matching_labels[:10]]

    @app.callback(
        Output(f'{prefix}-violin-plot2', 'figure'),
        [Input(f'{prefix}-violin2-gene-selection', 'value'),
         Input(f'{prefix}-multi-class-selection', 'value'),
         Input(f'{prefix}-binary-selection', 'value'),
         Input(f'{prefix}-p-value-method', 'value'),
         Input(f'{prefix}-show-box2', 'value'),
        Input(f'{prefix}-show-scatter2', 'value') ,
         Input(f'{prefix}-violin2-log-or-zscore', 'value'),
         Input(f'{prefix}-violin2-group-selection', 'value')]
    )
    def update_violin2(violin_gene_selection, multi_annotation, binary_annotation, p_value_method, 
                       show_box2, show_points, transformation, labels):
        p_value = None if p_value_method == 'none' else p_value_method
        return plot_violin2(
            adata,
            key=violin_gene_selection,
            groupby=multi_annotation,
            hue=binary_annotation,
            transformation=transformation,
            show_box='show' in show_box2 if show_box2 else False,
            show_points='show' in show_points if show_points else False,
            p_value=p_value,
            labels=labels
        )

    @app.callback(
        Output(f'{prefix}-dotplot', 'figure'),
        [Input(f'{prefix}-single-cell-genes-selection', 'value'),
         Input(f'{prefix}-single-cell-annotation-dropdown', 'value'),
         Input(f'{prefix}-single-cell-label-selection', 'value'),
         Input(f'{prefix}-plot-type-switch', 'value'),
         Input(f'{prefix}-aggregation-type', 'value'),
        Input(f'{prefix}-dotplot-log-or-zscore', 'value'),
         Input(f'{prefix}-dotplot-standardization', 'value'),
         Input(f'{prefix}-dotmatrix-color-map-dropdown', 'value')]
    )
    def update_dotplot(selected_genes, selected_annotation, selected_labels, plot_type, aggregation_type, 
                       transformation, standardization, color_map):
        filtered_adata = filter_data(adata, selected_annotation, selected_labels)
        return plot_dot_matrix(
            filtered_adata,
            selected_genes,
            selected_annotation,
            selected_labels,
            aggregation=aggregation_type,
            transformation=transformation,
            standardization=standardization,
            color_map=color_map,
            plot_type=plot_type
        )

    @app.callback(
        Output(f'{prefix}-stacked-bar-plot', 'figure'),
        [Input(f'{prefix}-x-meta', 'value'),
         Input(f'{prefix}-y-meta', 'value'),
         Input(f'{prefix}-norm-box', 'value'),
         Input(f'{prefix}-discrete-color-map-dropdown', 'value')]
    )
    def update_stacked_bar(x_meta, y_meta, norm, discrete_color_map):
        if x_meta == y_meta:
            raise PreventUpdate
        
        # Create color map if discrete color map is selected
        color_map = None
        if discrete_color_map:
            color_map = palette_json["color_palettes"][discrete_color_map]
        
        return plot_stacked_bar(x_meta, y_meta, norm, adata, color_map=color_map)
