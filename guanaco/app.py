import dash
import dash_bootstrap_components as dbc

# Initialize the app
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.LUX, "/assets/scientific_style.css"],
    suppress_callback_exceptions=True
) 