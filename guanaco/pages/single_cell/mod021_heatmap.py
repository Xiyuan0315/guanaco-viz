# from dash import dcc, html
# from config import common_config
# import dash_bootstrap_components as dbc
# import plotly.express as px
# import dash_draggable

# # Available color scales
# colorscales = px.colors.named_colorscales()

# def generate_heatmap_layout(adata, prefix):
#     """
#     Generate layout for the heatmap module.
#     """
    

#     heatmap_plot_draggable = html.Div([
#         dash_draggable.GridLayout(
#             id=f'{prefix}-heatmap-draggable',
#             className='grid-layout-no-border',
#             children=[
#                 html.Div(
#                     children=dcc.Loading(
#                         id=f"{prefix}-loading-heatmap",
#                         type="circle",
#                         children=[
#                             dcc.Graph(
#                                 id=f'{prefix}-heatmap',
#                                 config=common_config,
#                                 responsive=True,
#                                 style={
#                                     "min-height": "0",
#                                     "flex-grow": "1",
#                                     "border": "none",
#                                     "box-shadow": "none"
#                                 }
#                             )
#                         ],
#                         style={
#                             "height": "100%",
#                             "width": "100%",
#                             "display": "flex",
#                             "flex-direction": "column"
#                         }
#                     ),
#                     style={
#                         "height": '100%',
#                         "width": '100%',
#                         "display": "flex",
#                         "flex-direction": "column",
#                         "flex-grow": "0",
#                         "border": "none",
#                         "box-shadow": "none"
#                     }
#                 )
#             ],
#             isResizable=True,
#             isDraggable=True,
#             height=30,
#             gridCols=12,
#             style={
#                 "min-width": "300px",
#                 "min-height": "300px",
#                 "background": "transparent",
#                 "padding": "0px"
#             }
#         )
#     ], style={'height': '100%', 'display': 'flex', 'flex-direction': 'column', 'flex-grow': '1'})


#     heatmap_layout = html.Div([

# ,

#     html.Div(style={'marginTop': '20px'}),
#     html.Div([
#         heatmap_plot_draggable
#     ], style={'flex-grow': '1', 'min-height': '400px', 'display': 'flex', 'flex-direction': 'column'})
# ], style={'padding': '20px', 'display': 'flex', 'flex-direction': 'column', 'height': '100%'})


#     return heatmap_layout



from dash import dcc, html
import dash_bootstrap_components as dbc
import plotly.express as px
import dash_draggable
import json
from pathlib import Path
from guanaco.config import common_config

# Load color palettes
cvd_color_path = Path(__file__).parent / "cvd_color.json"
with open(cvd_color_path, "r") as f:
    palette_json = json.load(f)
palette_names = list(palette_json["color_palettes"].keys())


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
            value='None',
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

    # --- Secondary Annotation Color Map ---
    secondary_annotation_colormap_dropdown = html.Div([
        html.Label(
            'Secondary Annotation ColorMap:',
            style={'fontWeight': 'bold', 'marginBottom': '5px'}
        ),
        dcc.Dropdown(
            id=f'{prefix}-heatmap-secondary-colormap-dropdown',
            options=[{'label': name, 'value': name} for name in palette_names],
            value='tab20',  # Default to tab20 as you requested
            placeholder="Select colormap for secondary annotation",
            clearable=True,
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
    
    # dash_draggable.GridLayout(
    # id=f'{prefix}-draggable-heatmap',
    # className='grid-layout-no-border',
    # children=[
    #     html.Div([
    #         dcc.Graph(
    #             id=f'{prefix}-heatmap',
    #             config=common_config,
    #             responsive=True,
    #             style={
    #                 'min-height': '0',
    #                 'flex-grow': '1',
    #                 'backgroundColor': 'transparent'  # ensure the graph itself is transparent
    #             }
    #         )
    #     ],
    #         style={
    #             'backgroundColor': 'transparent', 
    #             'border': 'none', 
    #             'boxShadow': 'none',
    #             'height': '100%',
    #             'width': '100%',
    #             'display': 'flex',
    #             'flex-direction': 'column',
    #             'flex-grow': '0'  # ensure the div takes up available space
    #         })  # style of the div wrapper
    #     ],
    #     isResizable=True,
    #     isDraggable=True,
    #     height=30,
    #     gridCols=12,  # use gridCols for proper resizing
    #     style={
    #         'backgroundColor': 'transparent',
    #         'padding': '0px',
    #         'border': 'none',
    #         'boxShadow': 'none',
    #         "min-width": "300px",
    #         "min-height": "300px",
    #     }
    # )
    # --- Final layout container ---
    heatmap_layout = html.Div([
        heatmap_transformation_selection,
        heatmap_secondary_dropdown,
        dotmatrix_color_map_dropdown,
        secondary_annotation_colormap_dropdown,
        draggable_container
    ], style={'padding': '20px', 'marginBottom': '15px'})

    return heatmap_layout
