# ehymet 0.1.1

## Major Changes
- Implemented automatic variable selection functionality
  - Added new `vars_combinations = "auto"` option in `EHyClus`
  - Integrated `findCorrelation` and `select_var_ind` functions in `utils.R`
- Enhanced clustering validation metrics
  - Added Adjusted Rand Index (ARI) to `clustering_validation`
  - Added internal validations metrics such as Davies-Bouldin, Dunn and Silhouette
    indices available in `get_internal_clustering_criteria` function
  - Reordered parameters to `clusters` then `true_labels` for better intuitive use
  - Modified output to return results as a list

## Performance Improvements
- Added parallel execution support for the `generate_indices` function
- Implemented kmeans++ initialization for improved clustering
- Added support for additional parameters in `EHyClus`, `generate_indices`, and `funspline`

## Removals
- Removed the grid parameter (evaluation grid) as it wasn't affecting results
- Removed unused parameters for code cleanliness

## Documentation Updates
- Updated citations and added inst/CITATION
- Enhanced vignettes
  - Added section for automatic variable selection in second vignette
  - Updated examples and documentation

