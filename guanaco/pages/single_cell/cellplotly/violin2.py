import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from scipy.stats import ttest_ind, mannwhitneyu, f_oneway, kruskal
from scipy.stats import tukey_hsd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from guanaco.data_loader import color_config
from itertools import combinations
def calculate_p_values(df, method, groupby, hue, labels, comparison_mode='within'):
    """
    Calculate p-values for different comparison modes.
    
    comparison_mode: 'within' (compare hues within each group) or 'between' (compare groups within each hue)
    """
    if method == 'None' or method is None:
        return None

    p_values = {}
    unique_groups = df[groupby].unique()
    unique_hues = df[hue].unique()
    
    # For pairwise tests (2 levels)
    if method in ['mwu-test', 'ttest']:
        if comparison_mode == 'within':
            # Compare hues within each group
            for group in labels if labels else unique_groups:
                group_data = df[df[groupby] == group]
                if len(unique_hues) == 2:
                    # Pairwise comparison
                    group_1 = group_data[group_data[hue] == unique_hues[0]]['Expression']
                    group_2 = group_data[group_data[hue] == unique_hues[1]]['Expression']
                    if len(group_1) > 0 and len(group_2) > 0:
                        test_func = mannwhitneyu if method == 'mwu-test' else ttest_ind
                        kwargs = {'alternative': 'two-sided'} if method == 'mwu-test' else {'equal_var': False}
                        try:
                            _, p_val = test_func(group_1, group_2, **kwargs)
                            p_values[group] = p_val
                        except:
                            p_values[group] = np.nan
                else:
                    # Multiple comparisons - use tukey for pairwise
                    groups = []
                    labels_list = []
                    for h in unique_hues:
                        data = group_data[group_data[hue] == h]['Expression'].values
                        if len(data) > 0:
                            groups.append(data)
                            labels_list.append(h)
                    if len(groups) > 1:
                        p_values[group] = {'comparisons': [], 'p_values': []}
                        for (i, h1), (j, h2) in combinations(enumerate(labels_list), 2):
                            test_func = mannwhitneyu if method == 'mwu-test' else ttest_ind
                            kwargs = {'alternative': 'two-sided'} if method == 'mwu-test' else {'equal_var': False}
                            try:
                                _, p_val = test_func(groups[i], groups[j], **kwargs)
                                p_values[group]['comparisons'].append(f"{h1} vs {h2}")
                                p_values[group]['p_values'].append(p_val)
                            except:
                                pass
        else:  # between group comparison
            # Compare groups within each hue
            for hue_val in unique_hues:
                hue_data = df[df[hue] == hue_val]
                if len(unique_groups) == 2:
                    # Pairwise comparison
                    groups_list = sorted(unique_groups)
                    group_1 = hue_data[hue_data[groupby] == groups_list[0]]['Expression']
                    group_2 = hue_data[hue_data[groupby] == groups_list[1]]['Expression']
                    if len(group_1) > 0 and len(group_2) > 0:
                        test_func = mannwhitneyu if method == 'mwu-test' else ttest_ind
                        kwargs = {'alternative': 'two-sided'} if method == 'mwu-test' else {'equal_var': False}
                        try:
                            _, p_val = test_func(group_1, group_2, **kwargs)
                            p_values[hue_val] = p_val
                        except:
                            p_values[hue_val] = np.nan
                else:
                    # Multiple comparisons
                    groups = []
                    labels_list = []
                    for g in labels if labels else unique_groups:
                        data = hue_data[hue_data[groupby] == g]['Expression'].values
                        if len(data) > 0:
                            groups.append(data)
                            labels_list.append(g)
                    if len(groups) > 1:
                        p_values[hue_val] = {'comparisons': [], 'p_values': []}
                        for (i, g1), (j, g2) in combinations(enumerate(labels_list), 2):
                            test_func = mannwhitneyu if method == 'mwu-test' else ttest_ind
                            kwargs = {'alternative': 'two-sided'} if method == 'mwu-test' else {'equal_var': False}
                            try:
                                _, p_val = test_func(groups[i], groups[j], **kwargs)
                                p_values[hue_val]['comparisons'].append(f"{g1} vs {g2}")
                                p_values[hue_val]['p_values'].append(p_val)
                            except:
                                pass
    
    # For multi-group tests
    elif method in ['kw-test', 'anova']:
        if comparison_mode == 'within':
            # Compare multiple hues within each group - yields one p-value per group
            for group in labels if labels else unique_groups:
                group_data = df[df[groupby] == group]
                groups = []
                hue_names = []
                for h in unique_hues:
                    data = group_data[group_data[hue] == h]['Expression'].values
                    if len(data) > 0:
                        groups.append(data)
                        hue_names.append(h)
                if len(groups) > 1:
                    test_func = kruskal if method == 'kw-test' else f_oneway
                    try:
                        _, p_val = test_func(*groups)
                        p_values[str(group)] = p_val  # Ensure group is string
                    except Exception as e:
                        p_values[str(group)] = np.nan
        else:  # between group comparison
            # Compare multiple groups within each hue - yields one p-value per hue
            for hue_val in unique_hues:
                hue_data = df[df[hue] == hue_val]
                groups = []
                for g in labels if labels else unique_groups:
                    data = hue_data[hue_data[groupby] == g]['Expression'].values
                    if len(data) > 0:
                        groups.append(data)
                if len(groups) > 1:
                    test_func = kruskal if method == 'kw-test' else f_oneway
                    try:
                        _, p_val = test_func(*groups)
                        p_values[hue_val] = p_val
                    except:
                        p_values[hue_val] = np.nan

    elif method == 'two-way-anova':
        formula = f'Expression ~ C({groupby}) * C({hue})'
        try:
            model = ols(formula, data=df).fit()
            anova_table = sm.stats.anova_lm(model, typ=2)
            p_values['groupby'] = anova_table.loc[f'C({groupby})', 'PR(>F)']
            p_values['hue'] = anova_table.loc[f'C({hue})', 'PR(>F)']
            interaction_term = f'C({groupby}):C({hue})'
            p_values['interaction'] = anova_table.loc[interaction_term, 'PR(>F)']
        except:
            p_values = {'error': 'Failed to compute two-way ANOVA'}

    return p_values


def assign_colors(hue_levels, hue_color_map=None):
    if hue_color_map is None:
        hue_color_map = {}
    default_colors = px.colors.qualitative.Plotly
    for i, level in enumerate(sorted(hue_levels)):
        if level not in hue_color_map:
            hue_color_map[level] = color_config[i] if i < len(color_config) else default_colors[i % len(default_colors)]
    return hue_color_map


def add_p_value_annotations(fig, p_values, method, df, labels, comparison_mode='within', x_labels=None, group_levels=None, hue_levels=None, x_order=None, groupby=None, hue=None):
    """Add p-value annotations to the figure based on comparison mode."""
    if not p_values:
        return
    
    y_max = max(df['Expression']) * 1.05
    
    # Simple p-values (single value per group/hue)
    if all(isinstance(v, (int, float, np.number)) or (isinstance(v, float) and np.isnan(v)) for v in p_values.values()):
        if comparison_mode == 'within':
            # Annotations above each group
            # For grouped violins, we need to find the center position for each group
            if x_labels and group_levels and hue_levels and x_order:
                for group, p_val in p_values.items():
                    if (labels is None or group in labels) and not np.isnan(p_val):
                        # Find all x positions that belong to this group
                        group_positions = []
                        # In within mode, find positions based on x_order
                        for i, (g, h) in enumerate(x_order):
                            if str(g) == str(group):
                                group_positions.append(i)
                        
                        if group_positions:
                            # Place annotation at the center of the group
                            x_pos = np.mean(group_positions)
                            x_min = min(group_positions) - 0.3
                            x_max = max(group_positions) + 0.3
                            
                            # Add a bracket line to show which violins are being compared
                            fig.add_shape(
                                type="line",
                                x0=x_min, x1=x_max,
                                y0=y_max * 0.98, y1=y_max * 0.98,
                                line=dict(color="DarkSlateGrey", width=1.5),
                                xref="x", yref="y"
                            )
                            
                            fig.add_annotation(
                                x=x_pos,
                                y=y_max * 1.02,
                                text=f"<b>{p_val:.3f}</b>" if p_val < 0.05 else f"{p_val:.3f}",
                                showarrow=False,
                                font=dict(size=10, color='DarkSlateGrey'),
                                xref="x"
                            )
            else:
                # Fall back to original behavior for split violins
                for group, p_val in p_values.items():
                    if (labels is None or group in labels) and not np.isnan(p_val):
                        fig.add_annotation(
                            x=group,
                            y=y_max,
                            text=f"<b>{p_val:.3f}</b>" if p_val < 0.05 else f"{p_val:.3f}",
                            showarrow=False,
                            font=dict(size=10, color='DarkSlateGrey')
                        )
        else:  # between group comparison
            # Annotations above each hue (cluster)
            if x_labels and group_levels and hue_levels and x_order:
                for hue_val, p_val in p_values.items():
                    if not np.isnan(p_val):
                        # Find all x positions that belong to this hue
                        hue_positions = []
                        # In between mode, find positions based on x_order
                        for i, (g, h) in enumerate(x_order):
                            if str(h) == str(hue_val):
                                hue_positions.append(i)
                        
                        if hue_positions:
                            # Place annotation at the center of the hue group
                            x_pos = np.mean(hue_positions)
                            x_min = min(hue_positions) - 0.3
                            x_max = max(hue_positions) + 0.3
                            
                            # Add a bracket line to show which violins are being compared
                            fig.add_shape(
                                type="line",
                                x0=x_min, x1=x_max,
                                y0=y_max * 0.98, y1=y_max * 0.98,
                                line=dict(color="DarkSlateGrey", width=1.5),
                                xref="x", yref="y"
                            )
                            
                            fig.add_annotation(
                                x=x_pos,
                                y=y_max * 1.02,
                                text=f"<b>{p_val:.3f}</b>" if p_val < 0.05 else f"{p_val:.3f}",
                                showarrow=False,
                                font=dict(size=10, color='DarkSlateGrey'),
                                xref="x"
                            )
            else:
                # Create a summary annotation
                p_text = "\n".join([f"{k}: <b>{v:.3g}</b>" if v < 0.05 else f"{k}: {v:.3g}" 
                                    for k, v in p_values.items() if not np.isnan(v)])
                if p_text:
                    fig.add_shape(
                        type="path",
                        path="M 0.02 ,1.05 L 0.02,1.09 L 0.98,1.09 L 0.98,1.05",
                        xref="paper",
                        yref="paper",
                        line=dict(color="DarkSlateGrey", width=2)
                    )
                    fig.add_annotation(
                        text=p_text,
                        align="center",
                        showarrow=False,
                        xref="paper",
                        yref="paper",
                        x=0.5,
                        y=1.20,
                        font=dict(size=10, color="DarkSlateGrey"),
                        bgcolor="rgba(255, 255, 255, 0.8)"
                    )
    
    # Complex p-values (multiple comparisons)
    else:
        # Create a summary for multiple comparisons
        summary_text = []
        for key, value in p_values.items():
            if isinstance(value, dict) and 'comparisons' in value:
                summary_text.append(f"<b>{key}:</b>")
                for comp, p_val in zip(value['comparisons'], value['p_values']):
                    if p_val < 0.05:
                        summary_text.append(f"  {comp}: <b>{p_val:.3f}</b>")
                    else:
                        summary_text.append(f"  {comp}: {p_val:.3f}")
        
        if summary_text:
            # Place the text below the figure
            fig.add_annotation(
                text="<br>".join(summary_text),
                align="left",
                showarrow=False,
                xref="paper",
                yref="paper",
                x=0,
                y=-0.15,
                font=dict(size=9, color="DarkSlateGrey"),
                bgcolor="rgba(255, 255, 255, 0.8)",
                bordercolor="DarkSlateGrey",
                borderwidth=1
            )
    
    # Two-way ANOVA results
    if 'interaction' in p_values:
        anova_text = []
        if 'groupby' in p_values:
            anova_text.append(f"Main effect ({groupby}): {p_values['groupby']:.3f}")
        if 'hue' in p_values:
            anova_text.append(f"Main effect ({hue}): {p_values['hue']:.3f}")
        if 'interaction' in p_values:
            anova_text.append(f"Interaction: <b>{p_values['interaction']:.3f}</b>" 
                              if p_values['interaction'] < 0.05 
                              else f"Interaction: {p_values['interaction']:.3f}")
        
        fig.add_annotation(
            text="<br>".join(anova_text),
            align="left",
            showarrow=False,
            xref="paper",
            yref="paper",
            x=0.5,
            y=1.20,
            font=dict(size=10, color="DarkSlateGrey"),
            bgcolor="rgba(255, 255, 255, 0.8)"
        )


def plot_violin2(adata, key, groupby, hue, transformation,
                 show_box=True, show_points=True, p_value=None, labels=None,
                 hue_color_map=None, comparison_mode='within'):
    """
    Plot violin plots with support for multiple comparison modes.
    
    comparison_mode: 'within' (compare hues within each group) or 'between' (compare groups within each hue)
    """
    if not hue or not groupby or not key:
        return go.Figure()

    # Prepare data
    if key not in adata.var_names:
        return go.Figure()

    expression_data = adata[:, key].X.toarray().flatten()
    if transformation == 'log':
        expression_data = np.log1p(expression_data)
    elif transformation == 'z_score':
        expression_data = (expression_data - expression_data.mean()) / expression_data.std()

    df = pd.DataFrame({
        'Expression': expression_data,
        groupby: adata.obs[groupby],
        hue: adata.obs[hue],
        'Feature': key
    })

    hue_levels = sorted(df[hue].unique())
    group_levels = sorted(df[groupby].unique()) if labels is None else labels
    hue_color_map = assign_colors(hue_levels, hue_color_map)
    points_mode = 'all' if show_points else False

    fig = go.Figure()
    
    # Determine violin ordering based on comparison mode
    if comparison_mode == 'within':
        # Group-major ordering: all hues for each group stay together
        x_order = []
        x_labels = []
        for group in group_levels:
            for hue_val in hue_levels:
                x_order.append((group, hue_val))
                x_labels.append(f"{group}<br>{hue_val}")
    else:  # between
        # Hue-major ordering: all groups for each hue stay together
        x_order = []
        x_labels = []
        for hue_val in hue_levels:
            for group in group_levels:
                x_order.append((group, hue_val))
                x_labels.append(f"{group}<br>{hue_val}")
    
    # Create violins based on ordering
    if len(hue_levels) == 2 and comparison_mode == 'within':
        # Use split violins for within-group comparison with 2 hues
        for i, group in enumerate(group_levels):
            group_data = df[df[groupby] == group]
            for j, hue_val in enumerate(hue_levels):
                side = 'negative' if j == 0 else 'positive'
                sub_df = group_data[group_data[hue] == hue_val]
                if len(sub_df) > 0:
                    fig.add_trace(go.Violin(
                        x=[group] * len(sub_df),
                        y=sub_df['Expression'],
                        legendgroup=hue_val,
                        width=0.8,
                        scalemode='width',
                        scalegroup=str(group),
                        name=hue_val,
                        box_visible=show_box,
                        points=points_mode,
                        bandwidth=0.2,
                        spanmode='hard',
                        meanline_visible=True,
                        side=side,
                        line_color='DarkSlateGrey',
                        fillcolor=hue_color_map[hue_val],
                        hoveron='violins',
                        hoverinfo='y',
                        showlegend=(i == 0),
                        jitter=0.05,
                    ))
        fig.update_layout(violinmode='overlay')
    else:
        # Use grouped violins for all other cases
        for i, (group, hue_val) in enumerate(x_order):
            subset = df[(df[groupby] == group) & (df[hue] == hue_val)]
            if len(subset) > 0:
                # Show legend only for first occurrence of each hue
                show_legend = hue_val not in [t.legendgroup for t in fig.data]
                
                fig.add_trace(go.Violin(
                    x=[x_labels[i]] * len(subset),
                    y=subset['Expression'],
                    name=hue_val,
                    legendgroup=hue_val,
                    width=0.8,
                    meanline_visible=True,
                    points=points_mode,
                    box_visible=show_box,
                    bandwidth=0.2,
                    fillcolor=hue_color_map[hue_val],
                    line_color='DarkSlateGrey',
                    hoveron='violins',
                    hoverinfo='y',
                    spanmode='hard',
                    jitter=0.05,
                    showlegend=show_legend
                ))
        
        # Adjust layout for grouped violins
        if comparison_mode == 'within':
            fig.update_layout(violingap=0.3, violingroupgap=0.1)
        else:
            fig.update_layout(violingap=0.1, violingroupgap=0.3)

    # Calculate and add p-values
    if p_value and p_value != 'none':
        p_values = calculate_p_values(df, p_value, groupby, hue, labels=group_levels, comparison_mode=comparison_mode)
        # Pass x_labels only for grouped violins
        if len(hue_levels) == 2 and comparison_mode == 'within':
            add_p_value_annotations(fig, p_values, p_value, df, labels=group_levels, comparison_mode=comparison_mode,
                                    groupby=groupby, hue=hue)
        else:
            add_p_value_annotations(fig, p_values, p_value, df, labels=group_levels, 
                                    comparison_mode=comparison_mode, x_labels=x_labels, 
                                    group_levels=group_levels, hue_levels=hue_levels, x_order=x_order,
                                    groupby=groupby, hue=hue)

    # Check if we have complex p-values that will be displayed below
    has_complex_pvalues = p_value and p_value != 'none' and p_values and any(
        isinstance(v, dict) and 'comparisons' in v for v in p_values.values()
    )
    
    # Shared layout updates
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(size=12, color='DarkSlateGrey'),
        title=dict(
            text=f"Expression of {key} Across {groupby} by {hue}",
            x=0.5,
            xanchor='center',
            y=0.9,
            font=dict(size=16, color='DarkSlateGrey')
        ),
        yaxis_range=[0, max(df['Expression']) * 1.1],        
        xaxis_title=groupby if len(hue_levels) == 2 and comparison_mode == 'within' else "",
        yaxis_title='Expression',
        legend=dict(
            font=dict(size=10),
            title=dict(text=hue, font=dict(size=12))
        ),
        # Add bottom margin if we have complex p-values
        margin=dict(b=200 if has_complex_pvalues else 80)
    )

    if show_points:
        fig.update_traces(marker=dict(size=1.5), selector=dict(type='violin'))

    fig.update_yaxes(
        showgrid=False,
        showline=True, linewidth=2, linecolor='black',
        title=dict(text=f"{key}", font=dict(size=14, color='DarkSlateGrey'))
    )
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')

    return fig
