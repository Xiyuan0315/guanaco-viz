import json
import os
from pathlib import Path
import warnings
from dash import dcc, html, Input, Output, exceptions, State
from dash.exceptions import PreventUpdate
import dash
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import muon as mu
import anndata as ad

# Import visualization functions
from guanaco.pages.single_cell.cellplotly.embedding import plot_categorical_embedding, plot_continuous_embedding
from guanaco.pages.single_cell.cellplotly.heatmap1 import plot_heatmap1
from guanaco.pages.single_cell.cellplotly.heatmap2 import plot_heatmap2
from guanaco.pages.single_cell.cellplotly.violin1 import plot_violin1
from guanaco.pages.single_cell.cellplotly.violin2 import plot_violin2
from guanaco.pages.single_cell.cellplotly.stacked_bar import plot_stacked_bar
from guanaco.pages.single_cell.cellplotly.dotmatrix import plot_dot_matrix

# Import sub-layouts
from guanaco.pages.single_cell.mod021_heatmap import generate_heatmap_layout
from guanaco.pages.single_cell.mod022_violin import generate_violin_layout
from guanaco.pages.single_cell.mod023_dotplot import generate_dotplot_layout
from guanaco.pages.single_cell.mod024_stacked_bar import generate_stacked_bar_layout

# Import configs
from guanaco.config import scatter_config
from guanaco.data_loader import color_config

warnings.filterwarnings('ignore', message='.*observed=False.*')

# Load color palettes
cvd_color_path = Path(__file__).parent / "cvd_color.json"
with open(cvd_color_path, "r") as f:
    palette_json = json.load(f)
palette_names = list(palette_json["color_palettes"].keys())


# ============= Scatter Plot Functions =============

def initialize_scatter_components(adata):
    embedding_prefixes = {
        "X_umap": "UMAP", "X_pca": "PCA", "X_tsne": "t-SNE", "X_diffmap": "DiffMap",
        "X_phate": "PHATE", "X_draw_graph_fa": "FA"
    }
    
    obsm_list = list(adata.obsm.keys())
    embedding_columns = {
        key: [f'{embedding_prefixes.get(key, key.upper())}{i+1}' for i in range(adata.obsm[key].shape[1])]
        for key in obsm_list
    }
    default_embedding = obsm_list[-1]
    default_columns = embedding_columns[default_embedding]
    
    return obsm_list, embedding_columns, default_embedding, default_columns


def create_control_components(adata, prefix):
    obsm_list, _, default_embedding, default_columns = initialize_scatter_components(adata)
    clustering_dropdown = dcc.Dropdown(
        id=f'{prefix}-clustering-dropdown',
        options=[{'label': key.replace("X_", "").upper(), 'value': key} for key in obsm_list],
        value='X_umap' if 'X_umap' in obsm_list else default_embedding,
        placeholder="Select Clustering Method",
        style={'marginBottom': '15px','fontSize': '14px'},
        clearable=False
    )
    
    coordinates_dropdowns = html.Div([
        html.Label("X-axis:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        dcc.Dropdown(id=f'{prefix}-x-axis', options=[{'label': col, 'value': col} for col in default_columns], value=default_columns[0]),
        html.Label("Y-axis:", style={'fontWeight': 'bold', 'marginTop': '15px', 'marginBottom': '5px'}),
        dcc.Dropdown(id=f'{prefix}-y-axis', options=[{'label': col, 'value': col} for col in default_columns], value=default_columns[1])
    ], style={'marginBottom': '15px','fontSize': '14px'})
    
    return clustering_dropdown, coordinates_dropdowns


def generate_annotation_dropdown(anno_list, prefix):
    return dcc.Dropdown(id=f'{prefix}-annotation-dropdown', 
    options=[{'label': label, 'value': label} for label in anno_list],
    placeholder="Search and select an annotation...", 
    value = anno_list[len(anno_list)//2],
    style={'marginBottom': '15px'})


def generate_scatter_gene_selection(init_gene_list, prefix):
    return dcc.Dropdown(id=f'{prefix}-scatter-gene-selection', options=[{'label': label, 'value': label} for label in init_gene_list], 
    value = init_gene_list[0], placeholder="Search and select a gene...", style={'marginBottom': '15px'})


def scatter_layout(adata, prefix):
    scatter_transformation_selection = html.Div([
        dbc.RadioItems(
            id=f'{prefix}-scatter-log-or-zscore',
            options=[
                {'label': 'None', 'value':None},
                {'label': 'Log', 'value': 'log'},
                {'label': 'Z-score', 'value': 'z_score'},
            ],
            value='log',
            inline=True,
            style={'fontSize': '14px'} 
        )
    ])
    
    scatter_order_selection = html.Div([
        dbc.RadioItems(
            id=f'{prefix}-plot-order',
            options=[
                {'label': 'Max-1st', 'value': 'max'},
                {'label': 'Min-1st', 'value': 'min'},
                {'label': 'Original', 'value': 'original'},
                {'label': 'Random', 'value': 'random'},
            ],
            value='max',
            inline=True,
            style={'fontSize': '14px'}
        )
    ])
    
    colorscales = px.colors.named_colorscales()
    color_map_dropdown = dcc.Dropdown(
        id=f"{prefix}-scatter-color-map-dropdown",
        options=[{"label": scale, "value": scale} for scale in colorscales],
        value="viridis",
        style={"marginBottom": "10px"},
        clearable = False
    )
    
    color_map_discrete_dropdown = dcc.Dropdown(
        id=f"{prefix}-scatter-discrete-color-map-dropdown",
        options=[{"label": name, "value": name} for name in palette_names],
        value=None,
        placeholder="Default color",
        style={"marginBottom": "10px"},
        clearable=True,
        className="custom-dropdown",
    )
    
    marker_size_slider = dcc.Slider(
        id=f'{prefix}-marker-size-slider',
        min=1,
        max=10,
        value=5,
        marks = None,
        tooltip={"placement": "bottom", "always_visible": True},
        className="dbc-slider"
    )
    
    opacity_slider = dcc.Slider(
        id=f"{prefix}-opacity-slider",
        min=0.1,
        max=1.0,
        value=1,
        marks = None,
        tooltip={"placement": "bottom", "always_visible": True},
        className="dbc-slider"
    )
    
    scatter_legend_toggle = dbc.RadioItems(
        id=f'{prefix}-scatter-legend-toggle',
        options=[
            {'label': 'on data', 'value': 'on data'},
            {'label': 'right', 'value': 'right'}
        ],
        value='right',
        inline=True,
        style={'fontSize': '14px'} 
    )
    
    axis_toggle = dbc.RadioItems(
        id=f'{prefix}-axis-toggle',
        options=[
            {'label': 'Show', 'value': True},
            {'label': 'Hide', 'value': False}
        ],
        value=True,
        inline=True,
        style={'fontSize': '14px'}
    )
    
    graphic_control = html.Div(
        id=f"{prefix}-controls-container",
        children=[
            html.Div(
                [
                    html.Label("Plot order: ", className="control-label"),
                    scatter_order_selection,
                ],
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Continous ColorMap:",  className="control-label"),
                    color_map_dropdown,
                ],
                className="dbc",
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Discrete ColorMap:",  className="control-label"),
                    color_map_discrete_dropdown,
                ],
                className="dbc",
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Marker Size:",  className="control-label"),
                    marker_size_slider,
                ],
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Opacity:",  className="control-label"),
                    opacity_slider,
                ],
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Legend Location:",  className="control-label"),
                    scatter_legend_toggle,
                ],
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Axis:",  className="control-label"),
                    axis_toggle,
                ],
                style={'marginBottom': '15px'}
            ),
        ],
        style={'display': 'none'}  # hidden by default
    )
    
    unique_counts = adata.obs.nunique()
    unique_counts = unique_counts[unique_counts<100]
    unique_counts_sorted = unique_counts.sort_values()
    anno_list = unique_counts_sorted.index.tolist()
    
    clustering_dropdown, coordinates_dropdowns = create_control_components(adata, prefix)
    layout = html.Div([
        # Add Store component to hold selected cells data
        dcc.Store(id=f'{prefix}-selected-cells-store'),
        # Add a div to show selection status
        html.Div(id=f'{prefix}-selection-status', style={'textAlign': 'center', 'marginBottom': '10px'}),
        
        # Wrapper div to stabilize layout
        html.Div([
            dbc.Row([
            # Left column: Controls
            dbc.Col(
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("Dimension Reduction:", className = "control-label"),
                            clustering_dropdown,
                        ],
                        className="dbc",
                        style={'marginBottom': '15px'}
                    ),
                    html.Div(
                        [
                            html.Div(coordinates_dropdowns, id=f'{prefix}-coordinates-dropdowns')
                        ],
                        className="dbc",
                        style={'marginBottom': '15px'}
                    ),
                    html.Div(
                        [
                            html.Label("Transformation:",  className = "control-label"),
                            scatter_transformation_selection,
                        ],
                        style={'marginBottom': '15px'}
                    ),
                    html.Button("More controls", id=f"{prefix}-toggle-button", n_clicks=0, style={'marginBottom': '10px','border':'1px solid','borderRadius': '5px'}),
                    graphic_control,
                ],
            ),
            xs=12, sm=12, md=4, lg=4, xl=2, 
            style={"borderRight": "1px solid #ddd", "padding": "10px"}
        ),
        
        dbc.Col(
            html.Div(
                [
                    html.Label("Select Annotation:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                    generate_annotation_dropdown(anno_list=anno_list, prefix=prefix),
                    # Graph without loading component
                    dcc.Graph(
                        id=f'{prefix}-annotation-scatter', 
                        config=dict(scatter_config, **{'modeBarButtonsToAdd': ['select2d', 'lasso2d']}),
                        style={"height": "450px", "minHeight": "450px", "display": "block"}  # Match CSS min-height
                    ),
                ],
                className="dbc",
                style={'marginBottom': '20px'}
            ),
            xs=12, sm=12, md=4, lg=4, xl=5
        ),
        dbc.Col(
            html.Div(
                [
                    html.Label("Search Gene:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                    generate_scatter_gene_selection(init_gene_list=adata.var_names.to_list()[:10], prefix=prefix),
                    dcc.Graph(
                        id=f'{prefix}-gene-scatter', 
                        config=dict(scatter_config, **{'modeBarButtonsToAdd': ['select2d', 'lasso2d']}),
                        style={"height": "450px", "minHeight": "450px", "display": "block"}  # Match CSS min-height
                    ),
                ],
                className="dbc",
                style={'marginBottom': '20px','marginRight': '10px'}
            ),
            xs=12, sm=12, md=4, lg=4, xl=5
        ),
        ]),
        ], style={'position': 'relative'}),  # End wrapper div
        
        # Add a new row for the button - completely outside the scatter plot area
        dbc.Row([
            dbc.Col(
                # Empty column to align with controls
                width={"size": 2, "offset": 0}
            ),
            dbc.Col(
                html.Div([
                    dbc.Button(
                        "Update Other Plots with Selected Cells",
                        id=f"{prefix}-update-plots-button",
                        color="primary",
                        n_clicks=0,
                        size="md",
                        style={'width': '50%'}
                    ),
                    html.P(
                        "Select cells in the annotation scatter plot above, then click this button to filter all other plots",
                        style={'textAlign': 'center', 'marginTop': '10px', 'fontSize': '14px', 'color': '#666'}
                    )
                ], style={'marginTop': '20px', 'marginBottom': '20px'}),
                width={"size": 10}
            )
        ])
    ])
    
    return layout


# ============= Other Plots Functions =============

def filter_data(adata, annotation, selected_labels, selected_cells=None):
    """Filter data based on annotation labels and/or selected cells"""
    if selected_cells is not None and len(selected_cells) > 0:
        # Use selected cells if available
        adata_filtered = adata[selected_cells]
    elif selected_labels:
        # Otherwise use label filtering
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
    
    return html.Div([
        # Add a div to show cell filtering status
        html.Div(id=f'{prefix}-filter-status', style={'textAlign': 'center', 'marginBottom': '10px'}),
        dbc.Row([
            dbc.Col(
                generate_left_control(default_gene_markers, discrete_label_list, prefix),
                xs=12, sm=12, md=4, lg=4, xl=2
            ),
            dbc.Col(tabs, xs=12, sm=12, md=8, lg=8, xl=10)
        ], style={'marginBottom': '50px'})
    ])


# ============= Main Callback Functions =============

def single_cell_callbacks(app, adata, prefix):
    """Combined callback registration for both scatter and other plots"""
    
    # ===== Scatter Plot Callbacks =====
    
    @app.callback(
        Output(f"{prefix}-controls-container", "style"),
        Output(f"{prefix}-toggle-button", "children"),
        Input(f"{prefix}-toggle-button", "n_clicks"),
        prevent_initial_call=True
    )
    def toggle_controls(n_clicks):
        if n_clicks % 2 == 1:
            return {'display': 'block'}, "Hide controls"
        else:
            return {'display': 'none'}, "More controls"
    
    @app.callback(
    [Output(f'{prefix}-coordinates-dropdowns', 'children'),
    Output(f'{prefix}-x-axis', 'value'),
    Output(f'{prefix}-y-axis', 'value')],
    Input(f'{prefix}-clustering-dropdown', 'value')
    )
    def update_coordinates_dropdowns(selected_clustering):
        _, embedding_columns, _, _ = initialize_scatter_components(adata)
        selected_columns = embedding_columns[selected_clustering]
        options = [{'label': col, 'value': col} for col in selected_columns]
        x_value = selected_columns[0]
        y_value = selected_columns[1] if len(selected_columns) > 1 else selected_columns[0]
        
        return (
            html.Div([
                html.Div([
                    html.Label("X-axis:"),
                    dcc.Dropdown(id=f'{prefix}-x-axis', options=options, value=x_value, clearable=False,style={'fontSize': '14px'})
                ], style={'flex': '1', 'paddingRight': '10px'}),
                
                html.Div([
                    html.Label("Y-axis:"),
                    dcc.Dropdown(id=f'{prefix}-y-axis', options=options, value=y_value,clearable=False,style={'fontSize': '14px'})
                ], style={'flex': '1', 'paddingLeft': '10px'})
            ], style={'display': 'flex', 'marginBottom': '15px'}),
            x_value,
            y_value
        )
    
    @app.callback(
        Output(f'{prefix}-annotation-dropdown', 'options'),
        Input(f'{prefix}-annotation-dropdown', 'search_value')
    )
    def update_annotation_dropdown(search_value):
        if not search_value:
            raise exceptions.PreventUpdate
        label_list = adata.obs.keys().to_list()
        matching_labels = [label for label in label_list if search_value.lower() in label.lower()]
        return [{'label': label, 'value': label} for label in matching_labels[:10]]
    
    @app.callback(
        Output(f'{prefix}-scatter-gene-selection', 'options'),
        Input(f'{prefix}-scatter-gene-selection', 'search_value')
    )
    def update_scatter_gene_selection(search_value):
        if not search_value:
            raise exceptions.PreventUpdate
        gene_list = adata.var_names.to_list()
        matching_genes = [gene for gene in gene_list if search_value.lower() in gene.lower()]
        return [{'label': gene, 'value': gene} for gene in matching_genes[:20]]
    
    @app.callback(
        Output(f'{prefix}-annotation-scatter', 'figure'),
        [Input(f'{prefix}-clustering-dropdown', 'value'),
         Input(f'{prefix}-x-axis', 'value'),
         Input(f'{prefix}-y-axis', 'value'),
         Input(f'{prefix}-annotation-dropdown', 'value'),
        Input(f'{prefix}-scatter-gene-selection', 'value'),
         Input(f'{prefix}-marker-size-slider', 'value'),
         Input(f'{prefix}-opacity-slider', 'value'),
         Input(f'{prefix}-scatter-legend-toggle', 'value'),
         Input(f'{prefix}-axis-toggle', 'value'),
         Input(f'{prefix}-scatter-discrete-color-map-dropdown', 'value'),
         ]
    )
    def update_annotation_scatter(clustering_method, x_axis, y_axis,annotation, gene,marker_size, opacity, legend_show,axis_show,color_map):
        if not annotation:
            raise exceptions.PreventUpdate
        if color_map is None:
            color_map = color_config
        else:
            color_map = palette_json["color_palettes"][color_map]
        
        fig = plot_categorical_embedding(
            adata=adata,
            gene = gene,
            embedding_key=clustering_method,
            color=annotation,
            x_axis=x_axis,
            y_axis=y_axis,
            color_map=color_map,
            marker_size=marker_size,
            opacity=opacity,
            legend_show=legend_show,
            axis_show=axis_show,
        )
        
        # Enable selection mode and set height to match CSS
        fig.update_layout(
            dragmode='select',
            height=450,
            margin=dict(t=20, b=20, l=20, r=20)
        )
        
        return fig
    
    @app.callback(
        Output(f'{prefix}-gene-scatter', 'figure'),
        [Input(f'{prefix}-scatter-gene-selection', 'value'),
         Input(f'{prefix}-annotation-dropdown', 'value'),
         Input(f'{prefix}-clustering-dropdown', 'value'),
         Input(f'{prefix}-x-axis', 'value'),
         Input(f'{prefix}-y-axis', 'value'),
         Input(f'{prefix}-scatter-log-or-zscore', 'value'),
         Input(f'{prefix}-plot-order', 'value'),
         Input(f'{prefix}-scatter-color-map-dropdown', 'value'),
         Input(f'{prefix}-marker-size-slider', 'value'),
         Input(f'{prefix}-opacity-slider', 'value'),
         Input(f'{prefix}-annotation-scatter', 'relayoutData'),
         Input(f'{prefix}-axis-toggle', 'value'),
         ]
    )
    def update_gene_scatter(gene_name, annotation, clustering, x_axis, y_axis, transformation, order, color_map, marker_size, opacity, annotation_relayout, axis_show):
        if not gene_name:
            raise exceptions.PreventUpdate
        
        fig = plot_continuous_embedding(
            adata=adata,
            embedding_key=clustering,
            color=gene_name,
            x_axis=x_axis,
            y_axis=y_axis,
            transformation = transformation,
            order = order,
            color_map=color_map or 'Viridis',
            marker_size=marker_size,
            opacity=opacity,
            annotation = annotation,
            axis_show=axis_show,
        )
        
        # Set height to match CSS
        fig.update_layout(
            height=450,
            margin=dict(t=20, b=20, l=20, r=20)
        )
        
        # Apply zoom from annotation scatter plot
        if annotation_relayout and ('xaxis.range[0]' in annotation_relayout and 'yaxis.range[0]' in annotation_relayout):
            x_range = [annotation_relayout['xaxis.range[0]'], annotation_relayout['xaxis.range[1]']]
            y_range = [annotation_relayout['yaxis.range[0]'], annotation_relayout['yaxis.range[1]']]
            fig.update_layout(
                xaxis=dict(range=x_range), 
                yaxis=dict(range=y_range))
        
        return fig
    
    # ===== Cell Selection Callback =====
    @app.callback(
        Output(f'{prefix}-selected-cells-store', 'data'),
        Output(f'{prefix}-selection-status', 'children'),
        [Input(f'{prefix}-update-plots-button', 'n_clicks')],
        [State(f'{prefix}-annotation-scatter', 'selectedData'),
         State(f'{prefix}-annotation-dropdown', 'value')],
        prevent_initial_call=True
    )
    def store_selected_cells(n_clicks, selected_data, current_annotation):
        """Store indices of selected cells from annotation scatter plot when button is clicked"""
        if n_clicks == 0 or not selected_data or not selected_data.get('points'):
            return None, ""
        
        # Extract cell indices from selected points
        selected_points = selected_data['points']
        
        # Get the actual cell indices
        selected_indices = []
        for point in selected_points:
            # The curveNumber tells us which trace (category) the point belongs to
            curve_number = point.get('curveNumber', 0)
            point_number = point.get('pointNumber', 0)
            
            # Get all unique categories in order
            unique_categories = sorted(adata.obs[current_annotation].unique())
            selected_category = unique_categories[curve_number]
            
            # Get all cells in this category
            category_mask = adata.obs[current_annotation] == selected_category
            category_indices = adata.obs.index[category_mask].tolist()
            
            # The point_number is the index within this category
            if point_number < len(category_indices):
                selected_indices.append(category_indices[point_number])
        
        if selected_indices:
            n_selected = len(selected_indices)
            status_msg = dbc.Alert(
                f"✓ {n_selected} cells selected from {current_annotation}. Other plots updated.",
                color="success",
                dismissable=True,
                duration=4000
            )
            
            return selected_indices, status_msg
        else:
            return None, ""
    
    # ===== Other Plots Callbacks =====
    
    @app.callback(
        Output(f'{prefix}-violin2-group-selection', 'options'),
        Output(f'{prefix}-violin2-group-selection', 'value'),
        [Input(f'{prefix}-multi-class-selection', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data')]
    )
    def update_group_labels(selected_column, selected_cells):
        # Filter adata if cells are selected
        if selected_cells:
            filtered_adata = adata[selected_cells]
            unique_labels = sorted(filtered_adata.obs[selected_column].dropna().unique(), key=str)
        else:
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
        [Input(f'{prefix}-single-cell-annotation-dropdown', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data')]
    )
    def update_labels_based_on_annotation(selected_annotation, selected_cells):
        # Filter adata if cells are selected
        if selected_cells:
            filtered_adata = adata[selected_cells]
            unique_labels = filtered_adata.obs[selected_annotation].unique().tolist()
        else:
            unique_labels = adata.obs[selected_annotation].unique().tolist()
        
        label_options = [{'label': label, 'value': label} for label in sorted(unique_labels)]
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
         Input(f'{prefix}-discrete-color-map-dropdown', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data')]
    )
    def update_heatmap(selected_genes, selected_annotation, selected_labels, transformation, heatmap_color, secondary_annotation, discrete_color_map, selected_cells):
        filtered_adata = filter_data(adata, selected_annotation, selected_labels, selected_cells)
        
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
         Input(f'{prefix}-discrete-color-map-dropdown', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data')]
    )
    def update_violin1(selected_tab, selected_genes, selected_annotation, 
                       selected_labels, transformation, show_box_plot, show_scatter1, discrete_color_map, selected_cells):
        if selected_tab == 'violin-tab':
            # Use discrete color map if selected, otherwise use default
            color_map = None
            if discrete_color_map:
                discrete_palette = palette_json["color_palettes"][discrete_color_map]
                unique_labels = sorted(adata.obs[selected_annotation].unique())
                color_map = {
                    label: discrete_palette[i % len(discrete_palette)] for i, label in enumerate(unique_labels)
                }
            
            # Filter data if cells are selected
            filtered_adata = filter_data(adata, selected_annotation, selected_labels, selected_cells)
            
            fig = plot_violin1(
                filtered_adata,
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
         Input(f'{prefix}-dotmatrix-color-map-dropdown', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data')]
    )
    def update_dotplot(selected_genes, selected_annotation, selected_labels, plot_type, aggregation_type, 
                       transformation, standardization, color_map, selected_cells):
        filtered_adata = filter_data(adata, selected_annotation, selected_labels, selected_cells)
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
         Input(f'{prefix}-discrete-color-map-dropdown', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data')]
    )
    def update_stacked_bar(x_meta, y_meta, norm, discrete_color_map, selected_cells):
        # Check if only one metadata column is available
        unique_counts = adata.obs.nunique()
        discrete_columns = unique_counts[unique_counts < 100].index.tolist()
        
        if len(discrete_columns) < 2:
            # For single metadata, create a single stacked bar showing composition
            single_meta = discrete_columns[0]
            
            # Get counts for each category
            if selected_cells:
                filtered_adata = adata[selected_cells]
            else:
                filtered_adata = adata
                
            value_counts = filtered_adata.obs[single_meta].value_counts()
            
            # Create color mapping
            if discrete_color_map:
                discrete_palette = palette_json["color_palettes"][discrete_color_map]
                colors = [discrete_palette[i % len(discrete_palette)] for i in range(len(value_counts))]
            else:
                from guanaco.data_loader import color_config
                colors = [color_config[i % len(color_config)] for i in range(len(value_counts))]
            
            if norm == 'prop':
                # Convert to proportions
                values = (value_counts / value_counts.sum() * 100).values
                y_label = 'Percentage (%)'
                text_template = '%{y:.1f}%'
            else:
                values = value_counts.values
                y_label = 'Cell Count'
                text_template = '%{y}'
            
            # Create stacked bar chart with one bar
            fig = go.Figure()
            
            y_position = 0
            for i, (category, value) in enumerate(zip(value_counts.index, values)):
                fig.add_trace(go.Bar(
                    x=['All Cells'],
                    y=[value],
                    base=y_position,
                    name=str(category),
                    marker_color=colors[i],
                    hovertemplate=f'{category}<br>Count: {value_counts.iloc[i]}<br>%{{y:.1f}}%<extra></extra>' if norm == 'prop' else f'{category}<br>Count: %{{y}}<extra></extra>'
                ))
                y_position += value
            
            fig.update_layout(
                title=f'Composition by {single_meta}',
                barmode='stack',
                showlegend=True,
                plot_bgcolor='white',
                paper_bgcolor='white',
                xaxis=dict(
                    title='',
                    showgrid=False,
                    showline=True,
                    linewidth=2,
                    linecolor='black'
                ),
                yaxis=dict(
                    title=y_label,
                    showgrid=False,
                    showline=True,
                    linewidth=2,
                    linecolor='black'
                ),
                bargap=0.5,
                height=400
            )
            
            return fig
        
        if x_meta == y_meta or y_meta == 'none':
            # Create an empty figure with a message
            fig = go.Figure()
            fig.add_annotation(
                text="Please select different metadata columns for X-axis and color.",
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                showarrow=False,
                font=dict(size=14),
                xanchor="center",
                yanchor="middle"
            )
            fig.update_layout(
                plot_bgcolor='white',
                paper_bgcolor='white',
                xaxis=dict(visible=False),
                yaxis=dict(visible=False)
            )
            return fig
        
        # Filter data if cells are selected
        filtered_adata = adata[selected_cells] if selected_cells else adata
        
        # Create fixed color mapping based on ALL categories in the original dataset
        # This ensures consistent colors even when filtering
        all_y_categories = sorted(adata.obs[y_meta].unique())
        
        if discrete_color_map:
            discrete_palette = palette_json["color_palettes"][discrete_color_map]
            fixed_color_map = {
                cat: discrete_palette[i % len(discrete_palette)] 
                for i, cat in enumerate(all_y_categories)
            }
        else:
            from guanaco.data_loader import color_config
            fixed_color_map = {
                cat: color_config[i % len(color_config)] 
                for i, cat in enumerate(all_y_categories)
            }
        
        return plot_stacked_bar(x_meta, y_meta, norm, filtered_adata, color_map=fixed_color_map)
    
    # Add callback to update filter status
    @app.callback(
        Output(f'{prefix}-filter-status', 'children'),
        Input(f'{prefix}-selected-cells-store', 'data')
    )
    def update_filter_status(selected_cells):
        if selected_cells:
            n_selected = len(selected_cells)
            n_total = adata.n_obs
            return dbc.Alert(
                f"Showing {n_selected} of {n_total} cells ({n_selected/n_total*100:.1f}%) based on scatter plot selection",
                color="info",
                dismissable=False,
                style={'margin': '0'}
            )
        return None