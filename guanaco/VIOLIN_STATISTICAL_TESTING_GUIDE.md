# Statistical Testing Guide for Violin Plots

## Overview
This guide describes the statistical testing framework implemented in the Guanaco violin plot visualization tool. The system supports four analysis modes, each designed for different experimental designs and research questions.

## Analysis Modes

### Mode 1: Single Metadata Comparison
**Use Case**: Compare gene expression across groups in a single categorical variable (e.g., comparing expression across different cell types).

**Visualization**: Simple violin plots, one per group

**Statistical Tests**:
- **2 groups**: Mann-Whitney U test (non-parametric) or t-test (parametric)
- **>2 groups**: Kruskal-Wallis test (non-parametric) or ANOVA (parametric)

**Example**: Comparing CD8A expression across T cell subtypes (CD4+, CD8+, Tregs)

---

### Mode 2: Faceted Comparison
**Use Case**: Compare a secondary variable within subgroups of a primary variable (e.g., comparing treatment effect within each cell type).

**Visualization**: 
- **2 groups in meta2**: Split violins (side-by-side comparison)
- **>2 groups in meta2**: Grouped violins

**Statistical Tests**: Applied within each facet (meta1 group)
- **2 groups**: Mann-Whitney U test or t-test
- **>2 groups**: Kruskal-Wallis test or ANOVA

**Example**: Comparing treatment vs control (meta2) within each cell type (meta1)

**P-value Display**: One p-value above each facet group

---

### Mode 3: Linear Model with Confounder
**Use Case**: Test the effect of a primary variable while controlling for a confounding variable (e.g., testing treatment effect while controlling for batch effects).

**Model**: `expression ~ meta1 + meta2`

**Visualization**: Same as Mode 2 (grouped violins)

**Statistical Output**: 
- P-value for meta1 (main effect)
- P-value for meta2 (confounder effect)

**Example**: Testing drug response (meta1) while controlling for patient sex (meta2)

**P-value Display**: Single line showing both p-values (e.g., "Treatment: 0.023, Sex: 0.156")

---

### Mode 4: Mixed Effects Model
**Use Case**: Test the effect of a primary variable while accounting for random effects from grouped observations (e.g., multiple cells from the same patient).

**Model**: `expression ~ meta1 + (1|meta2)`
- meta1: Fixed effect (e.g., treatment condition)
- meta2: Random effect (e.g., patient ID)

**Visualization**: Same as Mode 2 (grouped violins)

**Statistical Approach**:
1. **Primary**: Linear mixed-effects model
2. **Fallback**: If convergence fails, uses pseudobulk approach:
   - Averages expression per meta2 group
   - Applies t-test or ANOVA on aggregated data

**Example**: Testing treatment effect (meta1) while accounting for patient-to-patient variability (meta2)

---

## Test Selection

### Automatic Selection
When "Auto" is selected, the tool automatically chooses the appropriate test based on:
1. **Mode type**
2. **Number of groups**
3. **Data characteristics**

### Manual Override
Users can manually select from available tests:
- Mann-Whitney U Test
- T-test
- Kruskal-Wallis Test
- ANOVA
- Linear Model
- Mixed Model

---

## Statistical Considerations

### Non-parametric vs Parametric Tests
- **Default**: Non-parametric tests (Mann-Whitney U, Kruskal-Wallis)
- **Rationale**: Single-cell data often violates normality assumptions
- **Alternative**: Parametric tests available for normally distributed data

### Multiple Testing Correction
- Currently, p-values are reported without correction
- For multiple comparisons, consider applying:
  - Bonferroni correction
  - Benjamini-Hochberg FDR
  - Adjust significance threshold accordingly

### Sample Size Considerations
- **Cell-level analysis**: Large n, but cells are not independent
- **Pseudobulk approach**: Reduces to biological replicates
- **Mixed models**: Account for hierarchical structure

---

## Interpretation Guidelines

### P-value Display
- **Bold**: p < 0.05 (statistically significant)
- **Regular**: p ≥ 0.05
- **Format**: 3 significant figures (e.g., 0.0234)

### Effect Size
While not displayed, consider:
- **Biological significance**: Is the difference meaningful?
- **Expression levels**: Check absolute expression values
- **Consistency**: Verify across cell types/conditions

### Common Pitfalls
1. **Pseudoreplication**: Treating cells as independent when they're from the same sample
2. **Batch effects**: Use Mode 3 to control for known confounders
3. **Zero inflation**: Common in single-cell data, may affect test validity

---

## Best Practices

### Choosing the Right Mode
1. **Single variable**: Use Mode 1
2. **Two variables, interested in within-group differences**: Use Mode 2
3. **Confounding variable to control**: Use Mode 3
4. **Hierarchical data structure**: Use Mode 4

### Data Preparation
- **Filter**: Remove low-quality cells
- **Normalize**: Ensure proper normalization
- **Transform**: Log transformation often helpful (available in UI)

### Validation
- **Visual inspection**: Check if statistical results match visual differences
- **Biological replicates**: Ensure sufficient samples per group
- **Technical replicates**: Consider averaging or mixed models

---

## Examples by Research Question

### "Is gene X differentially expressed between conditions?"
- **Mode 1**: If comparing conditions directly
- **Mode 3**: If need to control for batch/other factors

### "Does treatment affect gene X expression differently in each cell type?"
- **Mode 2**: Compare treatment within each cell type
- **Result**: Separate p-value for each cell type

### "Is gene X affected by treatment, accounting for patient variability?"
- **Mode 4**: Treatment as fixed effect, patient as random effect
- **Benefit**: More accurate p-values when patients vary

### "Do male and female patients respond differently to treatment?"
- **Mode 3**: Both sex and treatment as fixed effects
- **Result**: Independent p-values for each factor

---

## Technical Details

### Implementation
- **Statistical libraries**: SciPy, Statsmodels
- **Mixed models**: statsmodels.formula.api.mixedlm
- **Convergence**: Maximum likelihood estimation with fallback options

### Performance
- **Optimized for**: Large single-cell datasets
- **Caching**: Gene expression data cached for efficiency
- **Computation**: Statistical tests computed on-demand

---

## Future Enhancements
1. Multiple testing correction options
2. Effect size reporting
3. Power analysis tools
4. Post-hoc pairwise comparisons
5. Interaction term support for Mode 3