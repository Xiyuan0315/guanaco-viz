# Violin Plot Modes Update Summary

## Overview
Updated the violin plot functionality to support four analysis modes with appropriate statistical tests, replacing the previous "Comparison Type" approach.

## Files Modified

### 1. **violin2_new.py** (New file)
- Created new plotting function `plot_violin2_new()` that supports mode-based analysis
- Implemented automatic test selection based on mode and data characteristics
- Added support for:
  - **Mode 1**: One metadata comparison (t-test/MWU for 2 groups, ANOVA/KW for >2)
  - **Mode 2**: Faceted plots with within-facet comparisons
  - **Mode 3**: Linear models with confounders (expression ~ meta1 + meta2)
  - **Mode 4**: Mixed models with random effects (expression ~ meta1 + (1|meta2))
- Includes graceful fallback to pseudobulk for mixed models if convergence fails

### 2. **mod022_violin.py** 
- Replaced old dropdowns with new structure:
  - `meta1-selection`: Primary metadata
  - `meta2-selection`: Secondary metadata (optional, with "None" option)
  - `mode-selection`: Analysis mode (Mode 1-4)
  - `test-method-selection`: Statistical test (with "Auto" as default)
- Added mode explanation helper text div
- Updated layout to show new dropdown arrangement

### 3. **single_cell_plots.py**
- Added new callbacks:
  - `update_mode_explanation()`: Shows helpful text explaining each mode
  - `update_meta2_state()`: Disables meta2 for Mode 1
  - `update_test_methods()`: Filters available tests based on mode and data
- Updated main violin plot callback to use `plot_violin2_new()`
- Modified group selection callback to work with meta1

## Key Features

### Automatic Test Selection
When "Auto" is selected, the system automatically chooses:
- **Mode 1**: MWU (2 groups) or KW (>2 groups)
- **Mode 2**: MWU (2 groups) or KW (>2 groups) within each facet
- **Mode 3**: Linear model
- **Mode 4**: Mixed model with pseudobulk fallback

### Visualization Logic
- **Mode 1**: Simple violin plot grouped by meta1
- **Mode 2**: 
  - Split violins if meta2 has 2 levels
  - Grouped violins if meta2 has >2 levels
- **Mode 3 & 4**: Same visualization as Mode 2 but with model-based statistics

### P-value Display
- **Mode 1**: Single p-value above plot
- **Mode 2**: P-values for each facet
- **Mode 3 & 4**: Model summary below plot including:
  - Effect p-values
  - R² (for linear model)
  - Convergence warnings (for mixed model)
  - Fallback method notification

## Usage Examples

1. **Simple comparison across conditions**:
   - Meta1: "condition", Meta2: None, Mode: 1
   - Auto-selects t-test or ANOVA based on number of conditions

2. **Compare treatment within each cell type**:
   - Meta1: "cell_type", Meta2: "treatment", Mode: 2
   - Shows split violins with p-values for each cell type

3. **Control for batch effects**:
   - Meta1: "condition", Meta2: "batch", Mode: 3
   - Linear model adjusting for batch

4. **Account for patient variability**:
   - Meta1: "condition", Meta2: "patient", Mode: 4
   - Mixed model with patient as random effect

## Migration Notes
- The old "Comparison Type" (within/between) is replaced by the mode system
- Previous plots can be recreated:
  - "Within group comparison" → Mode 2
  - "Between group comparison" → Mode 2 (swap meta1/meta2)
- The system maintains backward compatibility through the original `plot_violin2()` function



  1. Clearer Mental Model

  Before: "Within group" vs "Between group" comparisons were confusing
  - Users had to think about directionality of comparisons
  - Unclear which variable was primary

  Now: Four modes that map to real research questions
  - Mode 1: "I have one variable"
  - Mode 2: "I want to see effects within subgroups"
  - Mode 3: "I need to control for confounders"
  - Mode 4: "My data has hierarchical structure"

  2. Automatic Intelligence

  - Auto mode selects the right test based on your data
  - No need to remember "use Mann-Whitney for 2 groups, Kruskal-Wallis for
  >2"
  >- System prevents invalid choices (e.g., can't select t-test for >2
  >groups)

  3. Better Terminology

  - "Meta1" and "Meta2" clearly indicate primary vs secondary variables
  - Mode names describe the analysis type, not the comparison direction
  - Helpful explanatory text appears when selecting each mode

  4. Simplified Output

  - P-values shown directly above relevant groups
  - Single-line summaries for complex models
  - No unnecessary technical details cluttering the visualization

  From the Mathematical/Statistical Perspective

  1. Proper Handling of Experimental Designs

  Confounding Variables (Mode 3)
  Old: Might compare Treatment vs Control, ignoring batch effects
  New: expression ~ treatment + batch (controls for batch)
  This gives unbiased treatment effect estimates by accounting for known
  confounders.

  2. Hierarchical Data Structure (Mode 4)

  The Problem: Single-cell data violates independence assumptions
  - 1000 cells from Patient A are not 1000 independent observations
  - Standard tests give falsely low p-values

  The Solution: Mixed effects models
  expression ~ treatment + (1|patient)
  - Accounts for patient-level clustering
  - More accurate p-values
  - Fallback to pseudobulk if model fails

  3. Appropriate Test Selection

  The system now considers:
  - Data distribution: Non-parametric tests as default (appropriate for
    zero-inflated single-cell data)
  - Sample structure: Different tests for paired vs unpaired designs
  - Effect types: Fixed effects vs random effects

  4. Statistical Power Optimization

  Mode 2 (Faceted): Tests within homogeneous subgroups
  - Higher power by reducing variance
  - Detects cell-type-specific effects

  Mode 3 (Linear Model): Increases power by explaining variance
  - Removing batch effect variance → cleaner signal
  - More sensitive to true biological differences

  5. Avoiding Common Pitfalls

  Pseudoreplication: Mode 4 explicitly handles this
  Simpson's Paradox: Mode 2 reveals within-group patterns that might
  reverse overall trends
  Confounding: Mode 3 provides unbiased estimates

  Real-World Example

  Research Question: "Does drug X affect gene expression?"

  Old Approach:
  - Simple t-test: drug vs control
  - Problem: Ignores patient effects, batch effects

  New Approach Options:
  1. Mode 1: If all samples are independent
  2. Mode 3: Drug + Batch (controls for batch effects)
  3. Mode 4: Drug + (1|Patient) (accounts for multiple cells per patient)

  Result: More accurate p-values that reflect true biological effects, not
  technical artifacts.

  Statistical Rigor Summary

  The new system provides:
  1. Unbiased estimates through proper model specification
  2. Correct Type I error rates by handling dependencies
  3. Increased power through variance decomposition
  4. Flexibility to match analysis to experimental design
  5. Robustness through automatic fallbacks

  This isn't just a UI improvement—it's a fundamental upgrade in
  statistical methodology that produces more reliable and interpretable
  results for single-cell analysis.