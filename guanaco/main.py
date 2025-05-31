from dash import dcc, html, Output, Input, MATCH
from guanaco.app import app
from guanaco.layout import navbar, tab_content, footprint, guanaco_footer, description_layout, anndata_layout, igv_layout
from guanaco.pages.browser.gene_browser import gene_browser_callbacks
from guanaco.pages.single_cell.mod01_scatter import scatter_callback
from guanaco.pages.single_cell.mod02_other_plots import maker_vis_callback
from guanaco.data_loader import initialize_data
import muon as mu
import anndata as ad

mu.set_options(pull_on_update=False)

# Load datasets - initialize_data will use environment variables set by CLI
datasets = initialize_data()

# Utility function
def get_discrete_labels(adata: ad.AnnData, *, max_unique: int = 50) -> list[str]:
    nunique = adata.obs.nunique()
    return nunique[nunique < max_unique].sort_values().index.tolist()

# App layout
app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    navbar(datasets),
    html.Div(id="tabs-content", style={"paddingTop": "70px"}),
    footprint,
    guanaco_footer
])

# Register callbacks for scatter and other plots for each dataset
for name, dataset in datasets.items():
    # Register AnnData callbacks if adata exists
    if dataset.adata is not None:
        if isinstance(dataset.adata, mu.MuData):
            for mod in dataset.adata.mod.keys():
                mod_adata = dataset.adata.mod[mod]
                prefix = f"{name}-{mod}"
                scatter_callback(app, mod_adata, prefix)
                maker_vis_callback(app, mod_adata, prefix)
        else:
            prefix = name
            scatter_callback(app, dataset.adata, prefix)
            maker_vis_callback(app, dataset.adata, prefix)

    # Register genome browser callbacks if genome tracks exist
    if dataset.genome_tracks is not None and dataset.ref_track is not None:
        gene_browser_callbacks(app, dataset.genome_tracks, dataset.ref_track, dataset.name)

# Update main content when dataset tab changes
@app.callback(
    Output("tabs-content", "children"),
    Input("tabs-dataset", "active_tab")
)
def update_tab_content(tab):
    dataset = datasets[tab]
    return tab_content(dataset, tab)

# Update description layout
@app.callback(
    Output({"type": "description-layout-div", "index": MATCH}, "children"),
    Input("tabs-dataset", "active_tab")
)
def update_description_layout(active_tab):
    dataset = datasets[active_tab]

    return description_layout(dataset)


# Update AnnData layout
@app.callback(
    Output({"type": "ann-layout-div", "index": MATCH}, "children"),
    Input({"type": "modality-tabs", "index": MATCH}, "active_tab"),
    Input("tabs-dataset", "active_tab")
)
def update_anndata_layout(selected_modality, active_tab):
    dataset = datasets[active_tab]
    if dataset.adata is None:
        return html.Div("No AnnData available for this dataset", style={"padding": "20px"})
    
    adata = dataset.adata.mod[selected_modality] if isinstance(dataset.adata, mu.MuData) else dataset.adata
    label_list = get_discrete_labels(adata)
    prefix = f"{active_tab}-{selected_modality}" if isinstance(dataset.adata, mu.MuData) else active_tab
    return anndata_layout(adata, dataset.gene_markers or [], label_list, prefix)

# Update IGV layout when dataset changes
@app.callback(
    Output({"type": "igv-layout-div", "index": MATCH}, "children"),
    Input("tabs-dataset", "active_tab")
)
def update_igv_layout(active_tab):
    dataset = datasets[active_tab]
    if dataset.genome_tracks is None:
        return html.Div("No genome browser data available for this dataset", style={"padding": "20px"})
    
    session_names = list(dataset.genome_tracks.keys())
    return igv_layout(session_names, prefix=active_tab)

server = app.server

if __name__ == "__main__":
    app.run_server(host='0.0.0.0', debug=True, port=4399)
