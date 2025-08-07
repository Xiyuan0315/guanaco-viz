from dash import dcc, html
import dash_bootstrap_components as dbc
import plotly.express as px
import dash_draggable
from guanaco.pages.single_cell.config import common_config


# Get the available color scales
colorscales = px.colors.named_colorscales()

def generate_heatmap_layout(adata, prefix):
    label_list = adata.obs.columns.to_list()

    # --- Transformation selection ---
    heatmap_transformation_selection =   html.Div([
        html.Label(
            'Transformation:',
            style={'fontWeight': 'bold', 'marginBottom': '5px'}
        ),
        dbc.RadioItems(
            id=f'{prefix}-heatmap-transformation',
            options=[
                {'label': 'None', 'value': 'None'},
                {'label': 'Log', 'value': 'log'},
                {'label': 'Z-score (across cell)', 'value': 'z_score'}
            ],
            value='log',
            inline=True,
            style={'marginBottom': '10px'}
        )
    ])

    # --- Secondary dropdown ---
    heatmap_secondary_dropdown =  html.Div([
        html.Label('Secondary Annotation:', style={'fontWeight': 'bold', 'marginTop': '10px'}),
        dcc.Dropdown(
            id=f'{prefix}-heatmap-label-dropdown',
            options=[{'label': 'None', 'value': 'None'}] + [{'label': label, 'value': label} for label in label_list],
            value='None',
            style={'width': '200px'}
        )
    ], style={'marginBottom': '10px'})

    # --- Color map selection ---
    dotmatrix_color_map_dropdown = html.Div([
        html.Label(
            'Color Map:',
            style={'fontWeight': 'bold', 'marginBottom': '5px'}
        ),
        dcc.Dropdown(
            id=f'{prefix}-heatmap-colorscale-dropdown',
            options=[
                {'label': scale, 'value': scale} for scale in colorscales
            ],
            value='viridis', 
            clearable=False,
            style={'width': '200px', 'marginBottom': '10px'}
        )
    ])

    # --- Draggable container with plot ---
    draggable_container = dash_draggable.GridLayout(
        id=f'{prefix}-draggable-heatmap',
        className='grid-layout-no-border',
        children=[
            html.Div(children=[
               dcc.Graph(
                  id=f'{prefix}-heatmap',
                  config=common_config,
                  responsive=True,
                  style={
                        "min-height":"0",
                        "flex-grow":"1"
                  })],
               style={
                  "height":'100%',
                  "width":'100%',
                  "display":"flex",
                  "flex-direction":"column",
                  "flex-grow":"0"
               }),
      ])

    heatmap_layout = html.Div([
        heatmap_transformation_selection,
        heatmap_secondary_dropdown,
        dotmatrix_color_map_dropdown,
        draggable_container
    ], style={'padding': '20px', 'marginBottom': '15px'})

    return heatmap_layout
