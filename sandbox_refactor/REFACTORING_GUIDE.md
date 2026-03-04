# Workflow Refactoring Guide

## Overview
Your main workflow has been refactored into 6 modular sub-workflows for better organization, maintainability, and reusability.

## Directory Structure
```
.
├── main_refactored.nf          # New main workflow using sub-workflows
├── subworkflows/               # New directory for sub-workflows
│   ├── diversity_metrics.nf
│   ├── frequency_analysis.nf
│   ├── merging_subsampling.nf
│   ├── imgt_aa_frequency.nf
│   ├── aa_properties_testing.nf
│   └── rarefied_sharing.nf
├── modules/                    # Your existing modules (unchanged)
└── data/                       # Your existing data (unchanged)
```

## Sub-workflows Created

### 1. **DIVERSITY_METRICS_WORKFLOW** (`subworkflows/diversity_metrics.nf`)
**Purpose**: Extract and combine diversity metrics for productive and non-productive sequences

**Inputs**:
- `productive_ch`: Channel of productive sequences (tuple: sample, file)
- `nonproductive_ch`: Channel of non-productive sequences (tuple: sample, file)

**Outputs**:
- `productive_combined`: Combined productive diversity metrics
- `nonproductive_combined`: Combined non-productive diversity metrics

**Processes included**:
- EXTRACT_DIVERSITY_METRICS_PROD
- EXTRACT_DIVERSITY_METRICS_NONPROD
- COMBINE_DIVERSITY_METRICS_PROD
- COMBINE_DIVERSITY_METRICS_NONPROD

---

### 2. **FREQUENCY_ANALYSIS_WORKFLOW** (`subworkflows/frequency_analysis.nf`)
**Purpose**: Compute and combine frequency metrics for productive and non-productive sequences

**Inputs**:
- `productive_ch`: Channel of productive sequences
- `nonproductive_ch`: Channel of non-productive sequences

**Outputs**:
- `productive_combined`: Combined productive frequencies
- `nonproductive_combined`: Combined non-productive frequencies

**Processes included**:
- COMPUTE_FREQS_PROD
- COMPUTE_FREQS_NONPROD
- COMBINE_FREQS_PROD
- COMBINE_FREQS_NONPROD

---

### 3. **MERGING_SUBSAMPLING_WORKFLOW** (`subworkflows/merging_subsampling.nf`)
**Purpose**: Merge samples, perform subsampling, compute frequencies, and calculate overlaps

**Inputs**:
- `productive_files_ch`: Channel of productive files
- `index_columns`: List of index columns for merging
- `count_column`: Name of count column
- `N`: Number of sequences per subsample
- `M`: Number of subsamples

**Outputs**:
- `merged_file`: Merged productive sequences
- `subsampled_files`: Collection of subsampled files
- `combined_freqs`: Combined frequencies from subsampled data
- `overlaps`: Overlap analysis results
- `overlaps_subsampled`: Overlap analysis on subsampled data

**Processes included**:
- MERGE_AMINO_VFAM
- SUBSAMPLE_FROM_MERGED
- COMPUTE_FREQS_SUBSAMPLED
- COMBINE_FREQS_SUBSAMPLED
- COMPUTE_OVERLAPS
- COMPUTE_OVERLAPS_SUBSAMPLED

---

### 4. **IMGT_AA_FREQUENCY_WORKFLOW** (`subworkflows/imgt_aa_frequency.nf`)
**Purpose**: Comprehensive IMGT amino acid frequency analysis with 4 variants:
- With length, no Vfam
- Without length, no Vfam
- With length, with Vfam
- Without length, with Vfam

**Inputs**:
- `subsampled_files_ch`: Channel of subsampled files
- `min_len`: Minimum CDR3 length
- `max_len`: Maximum CDR3 length
- `annotation_file`: Path to annotation file
- `min_non_index`: Minimum non-index value for filtering
- `min_value`: Minimum value for filtering

**Outputs**:
- `combined_subs_rows_full`: Combined rows (with length, no Vfam)
- `combined_subs_rows_WL`: Combined rows (without length, no Vfam)
- `combined_subs_rows_vfam`: Combined rows (with length, with Vfam)
- `combined_subs_rows_WL_vfam`: Combined rows (without length, with Vfam)
- `test_results_full`, `test_results_WL`, `test_results_vfam`, `test_results_WL_vfam`: Test results for each variant

**Processes included**:
- COMPUTE_IMGT_AA_FREQS_SUBS
- COMBINE_IMGT_AA_FREQS_MED (4 instances)
- PREP_COMB_IMGT_TEST (4 instances)
- RUN_IMGT_TEST (4 instances)

---

### 5. **AA_PROPERTIES_TESTING_WORKFLOW** (`subworkflows/aa_properties_testing.nf`)
**Purpose**: Test amino acid properties (VHS and Kidera factors) across all 4 variants

**Inputs**:
- `combined_subs_rows_full`: Combined rows (with length, no Vfam)
- `combined_subs_rows_WL`: Combined rows (without length, no Vfam)
- `combined_subs_rows_vfam`: Combined rows (with length, with Vfam)
- `combined_subs_rows_WL_vfam`: Combined rows (without length, with Vfam)
- `prop_VHS`: Path to VHS properties file
- `prop_kidera`: Path to Kidera factors file
- `annotation_file`: Path to annotation file

**Outputs**:
- 8 test result channels (VHS and Kidera for each of 4 variants)

**Processes included**:
- AA_PROPERTIES_TEST (8 instances: 4 for VHS, 4 for Kidera)

---

### 6. **RAREFIED_SHARING_WORKFLOW** (`subworkflows/rarefied_sharing.nf`)
**Purpose**: Perform rarefied sharing analysis with different clonotype definitions

**Inputs**:
- `merged_file_ch`: Channel of merged productive file
- `annotation_file`: Path to annotation file
- `K`: K parameter for rarefaction
- `M`: M parameter for rarefaction

**Outputs**:
- `sharing_full`: Sharing analysis (aminoAcid + vFamilyName + jGeneName)
- `sharing_vfam`: Sharing analysis (aminoAcid + vFamilyName)
- `sharing_aa`: Sharing analysis (aminoAcid only)

**Processes included**:
- RAREFIED_SHARING (3 instances with different column configurations)

---

## Main Workflow Changes

The new `main_refactored.nf` has been drastically simplified:

**Before**: ~200+ lines with complex channel operations scattered throughout
**After**: ~147 lines with clean sub-workflow calls

### Key improvements:
1. ✅ **Modularity**: Each logical unit is now a separate, reusable sub-workflow
2. ✅ **Readability**: Main workflow shows high-level flow without implementation details
3. ✅ **Maintainability**: Changes to a specific analysis only affect its sub-workflow
4. ✅ **Testability**: Each sub-workflow can be tested independently
5. ✅ **Reusability**: Sub-workflows can be used in other pipelines

## How to Use

### Option 1: Use the refactored main workflow
```bash
nextflow run main_refactored.nf
```

### Option 2: Import specific sub-workflows in your own scripts
```groovy
include { DIVERSITY_METRICS_WORKFLOW } from './subworkflows/diversity_metrics.nf'
include { IMGT_AA_FREQUENCY_WORKFLOW } from './subworkflows/imgt_aa_frequency.nf'

workflow {
    // Your custom workflow using just the sub-workflows you need
    DIVERSITY_METRICS_WORKFLOW(productive_ch, nonproductive_ch)
}
```

## Migration Steps

To switch from your original workflow to the refactored version:

1. **Backup your original main workflow**:
   ```bash
   cp main.nf main_original.nf
   ```

2. **Create the subworkflows directory**:
   ```bash
   mkdir -p subworkflows
   ```

3. **Copy the sub-workflow files** to `subworkflows/`:
   - diversity_metrics.nf
   - frequency_analysis.nf
   - merging_subsampling.nf
   - imgt_aa_frequency.nf
   - aa_properties_testing.nf
   - rarefied_sharing.nf

4. **Replace your main.nf** with `main_refactored.nf`:
   ```bash
   cp main_refactored.nf main.nf
   ```

5. **Test the refactored workflow**:
   ```bash
   nextflow run main.nf -profile test  # or your test profile
   ```

## Parameters

All parameters remain the same as your original workflow:

```groovy
params.input = './data/collated_info.csv'
params.outdir = './results'
params.seqs_to_remove = './data/D20210208D_1-overrepresented-sequences.txt'
params.N = 10000 
params.M = 25
params.index_columns = ["aminoAcid", "vFamilyName", "jGeneName"]
params.min_len = 8
params.max_len = 25
params.annotation_file = './data/collated_info.csv'
params.min_non_index = 5
params.min_value = 1
params.prop_VHS = 'assets/VHSE.tsv'
params.prop_kidera = 'assets/kidera_factors.tsv'
params.K = 200
```

## Benefits Summary

| Aspect | Before | After |
|--------|--------|-------|
| **Lines in main workflow** | ~200+ | ~147 |
| **Number of files** | 1 monolithic | 1 main + 6 sub-workflows |
| **Code reusability** | Low | High |
| **Testing complexity** | High (must test everything) | Low (test sub-workflows independently) |
| **Readability** | Complex with nested logic | Clean, high-level flow |
| **Maintenance** | Difficult (find code in large file) | Easy (each sub-workflow is focused) |

## Troubleshooting

### Issue: Module not found errors
**Solution**: Ensure your `modules/` directory structure matches the include paths. All sub-workflows expect modules to be in `../modules/` relative to the `subworkflows/` directory.

### Issue: Parameter not defined
**Solution**: Check that all required parameters are defined in your main workflow or config file.

### Issue: Channel type mismatch
**Solution**: Verify that the output channels from one sub-workflow match the expected input format of the next sub-workflow in the chain.

## Next Steps

1. Review each sub-workflow file to understand its structure
2. Test the refactored workflow with your data
3. Consider adding documentation comments to each sub-workflow
4. Create unit tests for individual sub-workflows
5. Add process-level error handling where appropriate

---

**Created by**: Seqera AI
**Date**: 2026-03-01
**Nextflow DSL Version**: DSL2
