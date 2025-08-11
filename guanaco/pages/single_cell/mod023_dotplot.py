from dash import dcc, html
import dash_bootstrap_components as dbc
import plotly.express as px
from guanaco.config import common_config
import dash_draggable

# Get the available color scales
colorscales = px.colors.named_colorscales()

def generate_dotplot_layout(prefix):
    """
    Generate layout for dotplot/matrixplot visualization.

    Arguments:
        prefix: Unique ID prefix for layout components.
    """
    # --- Plot type switch ---
    plot_type_switch = html.Div([
        html.Label(
            'Plot Type:',
            style={'fontWeight': 'bold', 'marginBottom': '5px'}
        ),
        dbc.RadioItems(
            id=f'{prefix}-plot-type-switch',
            options=[
                {'label': 'Dotplot', 'value': 'dotplot'},
                {'label': 'Matrixplot', 'value': 'matrixplot'}
            ],
            value='dotplot',
            inline=True,
            style={'marginBottom': '15px'}
        )
    ])


    # --- Transformation selection ---
    dotplot_transformation_selection = html.Div([
        html.Label(
            'Transformation:',
            style={'fontWeight': 'bold', 'marginBottom': '5px'}
        ),
        dbc.RadioItems(
            id=f'{prefix}-dotplot-log-or-zscore',
            options=[
                {'label': 'None', 'value': 'None'},
                {'label': 'Log', 'value': 'log'}
            ],
            value='None',
            inline=True,
            style={'marginBottom': '10px'}
        )
    ])

    # --- Standardization selection ---
    dotplot_standardization_selection = html.Div([
        html.Label(
            'Standardization:',
            style={'fontWeight': 'bold', 'marginBottom': '5px'}
        ),
        dbc.RadioItems(
            id=f'{prefix}-dotplot-standardization',
            options=[
                {'label': 'None', 'value': 'None'},
                {'label': 'By variable', 'value': 'var'},
                {'label': 'By group', 'value': 'group'}
            ],
            value='None',
            inline=True,
            style={'marginBottom': '10px'}
        )
    ])

    # --- Color map selection ---
    dotmatrix_color_map_dropdown = html.Div([
        html.Label(
            'Color Map:',
            style={'fontWeight': 'bold', 'marginBottom': '5px'}
        ),
        dcc.Dropdown(
            id=f'{prefix}-dotmatrix-color-map-dropdown',
            options=[
                {'label': scale, 'value': scale} for scale in colorscales
            ],
            value='viridis',  # Plotly scales are usually capitalized (like Viridis, Plasma, etc.)
            clearable=False,
            style={'width': '200px', 'marginBottom': '10px'}
        )
    ])

    # --- Draggable container with plot ---
    draggable_container = dash_draggable.GridLayout(
    id=f'{prefix}-draggable-dotplot',
    className='grid-layout-no-border',
    children=[
        html.Div([
            dcc.Graph(
                id=f'{prefix}-dotplot',
                config=common_config,
                responsive=True,
                style={
                    'min-height': '0',
                    'flex-grow': '1',
                    'backgroundColor': 'transparent'  # ensure the graph itself is transparent
                }
            )
        ],
            style={
                'backgroundColor': 'transparent', 
                'border': 'none', 
                'boxShadow': 'none',
                'height': '100%',
                'width': '100%',
                'display': 'flex',
                'flex-direction': 'column',
                'flex-grow': '0'  # ensure the div takes up available space
            })  # style of the div wrapper
        ],
        isResizable=True,
        isDraggable=True,
        height=30,
        gridCols=12,  # use gridCols for proper resizing
        style={
            'backgroundColor': 'transparent',
            'padding': '0px',
            'border': 'none',
            'boxShadow': 'none'
        }
    )
    # --- Final layout container ---
    dotplot_layout = html.Div([
        dotplot_transformation_selection,
        dotplot_standardization_selection,
        dotmatrix_color_map_dropdown,
        plot_type_switch,
        draggable_container
    ], style={'padding': '20px', 'marginBottom': '15px'})

    return dotplot_layout
