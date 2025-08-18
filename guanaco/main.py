from dash import dcc, html, Output, Input, MATCH, State
from guanaco.app import app
from guanaco.layout import (
    navbar, tab_content, footprint, guanaco_footer, description_layout,
    anndata_layout, igv_layout, resize_tip_toast  # make sure this is imported
)
from guanaco.pages.browser.gene_browser import gene_browser_callbacks
from guanaco.pages.single_cell.single_cell_plots import single_cell_callbacks
from guanaco.data_loader import datasets  # Already loaded at module level in data_loader.py
import muon as mu
import anndata as ad

mu.set_options(pull_on_update=False)

# Note: datasets is already loaded when data_loader.py was imported
# No need to call initialize_data() again!

# Utility function
def get_discrete_labels(adata: ad.AnnData, *, max_unique: int = 50) -> list[str]:
    nunique = adata.obs.nunique()
    return nunique[nunique < max_unique].sort_values().index.tolist()

# App layout
# app.layout = html.Div([
#     dcc.Location(id='url', refresh=False),
#     navbar(datasets),
#     html.Div(id="tabs-content", style={"paddingTop": "70px"}),
#     footprint,
#     guanaco_footer
# ])
app.layout = html.Div([
    dcc.Location(id="url", refresh=False),
    dcc.Store(id="tip-store", storage_type="session", data={"shown": False}),
    navbar(datasets),
    resize_tip_toast(),  # 👈 add this line
    html.Div(id="tabs-content", style={"paddingTop": "70px"}),
    footprint,
    guanaco_footer,
])


# Register callbacks for scatter and other plots for each dataset
for name, dataset in datasets.items():
    # Register AnnData callbacks if adata exists
    if dataset.adata is not None:
        if isinstance(dataset.adata, mu.MuData):
            for mod in dataset.adata.mod.keys():
                mod_adata = dataset.adata.mod[mod]
                prefix = f"{name}-{mod}"
                single_cell_callbacks(app, mod_adata, prefix)
        else:
            prefix = name
            single_cell_callbacks(app, dataset.adata, prefix)

    # Register genome browser callbacks if genome tracks exist
    if dataset.genome_tracks is not None and dataset.ref_track is not None:
        gene_browser_callbacks(app, dataset.genome_tracks, dataset.ref_track, dataset.title)

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
    
    # Use JSON markers for RNA, automatic markers for other modalities
    if selected_modality == 'rna' and dataset.gene_markers is not None:
        # Use markers from JSON config for RNA
        modality_markers = dataset.gene_markers
    else:
        # Use first 6 variables from the modality-specific adata
        modality_markers = adata.var_names[:6].tolist() if adata else []
    
    return anndata_layout(adata, modality_markers, label_list, prefix)

@app.callback(
    [Output("tip-modal", "is_open"), Output("tip-store", "data")],
    [Input("url", "pathname"), Input("close-tip", "n_clicks")],
    [State("tip-modal", "is_open"), State("tip-store", "data")],
)
def toggle_tip(pathname, n_clicks, is_open, store):
    if store is None:
        store = {"shown": False}
    
    if n_clicks:
        return False, {"shown": True}
    
    if not store.get("shown", False):
        return True, {"shown": True}

    return False, store


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
