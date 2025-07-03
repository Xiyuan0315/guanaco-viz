import json
import os
from dash import dcc, html, Input, Output, exceptions, State
import dash_bootstrap_components as dbc
import plotly.express as px
from guanaco.pages.single_cell.cellplotly.embedding import plot_categorical_embedding, plot_continuous_embedding
from guanaco.config import scatter_config
from guanaco.data_loader import color_config



# Load cvd_color.json from the same directory as this module
cvd_color_path = os.path.join(os.path.dirname(__file__), "cvd_color.json")
with open(cvd_color_path, "r") as f:
    palette_json = json.load(f)
palette_names = list(palette_json["color_palettes"].keys())

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

### 2. Initialize Scatter Components ###

def scatter_layout(adata,prefix):



    def create_control_components(adata):

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

    def generate_annotation_dropdown(anno_list):

        """Initial annotation dropdown with limited options."""
        return dcc.Dropdown(id=f'{prefix}-annotation-dropdown', 
        options=[{'label': label, 'value': label} for label in anno_list],
        placeholder="Search and select an annotation...", 
        value = anno_list[len(anno_list)//2],
        style={'marginBottom': '15px'})

    def generate_scatter_gene_selection(init_gene_list):
        """Initial gene selection dropdown with limited options."""
        return dcc.Dropdown(id=f'{prefix}-scatter-gene-selection', options=[{'label': label, 'value': label} for label in init_gene_list], 
        value = init_gene_list[0], placeholder="Search and select a gene...", style={'marginBottom': '15px'})

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
                {'label': 'Max-1st', 'value': 'max'},  # None option with 'False' value
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
    # Dropdown for color map selection
    color_map_dropdown = dcc.Dropdown(
        id=f"{prefix}-scatter-color-map-dropdown",
        options=[{"label": scale, "value": scale} for scale in colorscales],
        value="viridis",  # Default colorscale
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


    # Slider for marker size
    marker_size_slider = dcc.Slider(
        id=f'{prefix}-marker-size-slider',
        min=1,
        max=10,
        value=5,  # Default marker size
        marks = None,
        tooltip={"placement": "bottom", "always_visible": True},
        className="dbc-slider"
    )

    # Slider for opacity
    opacity_slider = dcc.Slider(
        id=f"{prefix}-opacity-slider",
        min=0.1,
        max=1.0,
        value=1,  # Default opacity
        marks = None,
        tooltip={"placement": "bottom", "always_visible": True},
        className="dbc-slider"
    )
    # toggle for scatter plot legend
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
    # Toggle for axis show
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

    clustering_dropdown, coordinates_dropdowns = create_control_components(adata)
    layout = html.Div([
        # Add Store component to hold selected cells data
        dcc.Store(id=f'{prefix}-selected-cells-store'),
        
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
    # First annotation sencond gene scatter plot

    dbc.Col(
        html.Div(
            [
                html.Label("Select Annotation:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                generate_annotation_dropdown(anno_list=anno_list),
                dcc.Loading(
                    id=f"{prefix}-loading-annotaion-scatter",
                    type="circle",
                    children=dcc.Graph(id=f'{prefix}-annotation-scatter', config=scatter_config),
                    style={"height": "100%"},
                ),
                # Add button below the scatter plot
                html.Div([
                    dbc.Button(
                        "Update Other Plots with Selected Cells",
                        id=f"{prefix}-update-plots-button",
                        color="primary",
                        n_clicks=0,
                        style={'width': '100%', 'marginTop': '10px'}
                    ),
                    html.Div(id=f"{prefix}-selection-status", style={'textAlign': 'center', 'marginTop': '5px'})
                ]),
            ],
            className="dbc",
            style={'marginBottom': '20px'}
        ),
        xs=12, sm=12, md=4, lg=4, xl=5 # Full width on small screens, half on larger screens
    ),
    dbc.Col(
        html.Div(
            [
                html.Label("Search Gene:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                generate_scatter_gene_selection(init_gene_list=adata.var_names.to_list()[:10]),
                dcc.Loading(
                    id=f"{prefix}-loading-gene-scatter",
                    type="circle",
                    children=dcc.Graph(id=f'{prefix}-gene-scatter', config=scatter_config),
                    style={"height": "100%"},
                ),
            ],
            className="dbc",
            style={'marginBottom': '20px','marginRight': '10px'}
        ),
        xs=12, sm=12, md=4, lg=4, xl=5  # Full width on small screens, half on larger screens
    ),
        ])
    ])

    return layout



## 3. Dynamic Dropdown Callbacks and Scatter Plot Callbacks ###
def scatter_callback(app, adata,prefix):

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
        """Dynamically update annotation options based on search input."""
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
         Input(f'{prefix}-gene-scatter', 'relayoutData'),
         Input(f'{prefix}-scatter-legend-toggle', 'value'),
         Input(f'{prefix}-axis-toggle', 'value'),
         Input(f'{prefix}-scatter-discrete-color-map-dropdown', 'value'),
         ]
    )
    def update_annotation_scatter(clustering_method, x_axis, y_axis,annotation, gene,marker_size, opacity, gene_relayout,legend_show,axis_show,color_map):

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
        # Enable selection tools
        fig.update_layout(
            dragmode='select',
            selectdirection='diagonal'
        )
        if gene_relayout and ('xaxis.range[0]' in gene_relayout and 'yaxis.range[0]' in gene_relayout):
            x_range = [gene_relayout['xaxis.range[0]'], gene_relayout['xaxis.range[1]']]
            y_range = [gene_relayout['yaxis.range[0]'], gene_relayout['yaxis.range[1]']]
            fig.update_layout(
                xaxis=dict(range=x_range), 
                yaxis=dict(range=y_range))


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
    def update_gene_scatter(gene_name, annotation, clustering, x_axis, y_axis, transformation, order, color_map, marker_size, opacity, gene_relayout,axis_show):
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
        if gene_relayout and ('xaxis.range[0]' in gene_relayout and 'yaxis.range[0]' in gene_relayout):
            x_range = [gene_relayout['xaxis.range[0]'], gene_relayout['xaxis.range[1]']]
            y_range = [gene_relayout['yaxis.range[0]'], gene_relayout['yaxis.range[1]']]
            fig.update_layout(
                xaxis=dict(range=x_range), 
                yaxis=dict(range=y_range))


        return fig
    
    # Add callback to handle selected cells and button click
    @app.callback(
        [Output(f'{prefix}-selected-cells-store', 'data'),
         Output(f'{prefix}-selection-status', 'children')],
        [Input(f'{prefix}-update-plots-button', 'n_clicks')],
        [State(f'{prefix}-annotation-scatter', 'selectedData'),
         State(f'{prefix}-clustering-dropdown', 'value'),
         State(f'{prefix}-x-axis', 'value'),
         State(f'{prefix}-y-axis', 'value'),
         State(f'{prefix}-annotation-dropdown', 'value')]
    )
    def update_selected_cells(n_clicks, selected_data, clustering_method, x_axis, y_axis, annotation):
        if n_clicks == 0 or not selected_data or not selected_data.get('points'):
            return None, ""
        
        # Get the selected points
        selected_points = selected_data['points']
        
        # Extract cell indices from the selected points
        # The points contain the indices of the cells in the original adata
        selected_indices = []
        for point in selected_points:
            # The curveNumber tells us which trace (category) the point belongs to
            curve_number = point.get('curveNumber', 0)
            point_number = point.get('pointNumber', 0)
            
            # We need to reconstruct the actual cell index
            # Get all unique categories in order
            unique_categories = sorted(adata.obs[annotation].unique())
            selected_category = unique_categories[curve_number]
            
            # Get all cells in this category
            category_mask = adata.obs[annotation] == selected_category
            category_indices = adata.obs.index[category_mask].tolist()
            
            # The point_number is the index within this category
            if point_number < len(category_indices):
                selected_indices.append(category_indices[point_number])
        
        if selected_indices:
            n_selected = len(selected_indices)
            status_msg = html.Div([
                html.Span(f"✓ {n_selected} cells selected from ", style={'color': 'green'}),
                html.Span(f"{annotation}", style={'fontWeight': 'bold'}),
                html.Span(". Other plots will show only these cells.", style={'color': 'green'})
            ])
            return selected_indices, status_msg
        else:
            return None, ""