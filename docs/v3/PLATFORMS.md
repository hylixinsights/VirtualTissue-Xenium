# Platform inputs and units

All adapters return annotated units, exact feature symbols and a `measurement` record. Existing object labels are accepted with their declared source; independent manual annotation review is not required. Map those labels explicitly to model lineages. The prepared bundle fingerprints the selected counts, coordinates, declarations, ontology and signatures. Native Atera/Visium/Stereo readers have format-fixture tests; a public real-data benchmark on each platform remains required before claiming cross-platform validity.

| Platform | Native input implemented | Required declarations |
|---|---|---|
| Xenium | Existing 10x H5 + cells table + existing-label CSV | Organism checked against local SOURCE.json when present |
| Atera | `cell_feature_matrix.zarr.zip`, cell-by-feature CSR | Supplied barcode labels and lineage mapping; actual gene feature category if different from `Gene Expression` |
| Visium HD | `filtered_feature_bc_matrix.h5`, `spatial/tissue_positions.parquet`, scale JSON | Micrometres per pixel and bin size; lattice for unsegmented bins |
| Stereo-seq | GEM/TSV, optionally gzip; segmented CellID or pre-binned units | Explicit coordinate scale, unit type, complete gene universe or explicit observed-union limitation; segmented-cell centroids supplied in annotation CSV |
| scRNA/snRNA | Annotated H5AD or JSON | Well-mixed only; nucleus unit for snRNA |

Atera's actual output is Zarr, superseding the original plan's anticipated HDF5 format. The reader uses CSR row slices and per-cell centroid coordinates in micrometres, without importing vendor cluster/annotation predictions as reviewed lineages. [10x matrix specification](https://www.10xgenomics.com/support/software/atera-onboard-analysis/latest/analysis/output-zarr-cfmatrix).

Visium HD pixel coordinates are converted using a declared calibrated scale; no fixed conversion is guessed. The reader preserves measured bins and does not create a segmentation mask. Existing externally segmented objects can be provided through H5AD/JSON with a declared segmentation method. [10x spatial output specification](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/spatial-outputs).

Stereo-seq GEM does not provide a complete list of zero-count assayed genes. Supply `measurement.assayed_genes` for a panel-consistent score; exploratory use of `feature_universe: observed_GEM_union` is explicitly flagged. `geneName` is preferred when present, otherwise `geneID` is retained exactly. Ensembl-to-symbol conversion is never inferred. GEF requires export using vendor tools; a native GEF reader is not implemented. [STOmics file-format manual](https://cdn-newfile.stomics.tech/Stereo-seq_Analysis_Workflow_File_Format_A4.2.pdf).

## Example import

If an annotation table is needed, create a blank sheet from measured geometry and join the object's existing labels by exact ID:

```bash
vtx2 platform-template visium_hd path/to/square_008um annotations.csv \
  --measurement measurement.json --center 1200 1600 --max-cells 1000
```

The same command supports `atera` and `stereo_seq`. It writes measured IDs/XY and empty `cell_type`, `lineage` and `purity` columns; vendor predictions are not silently treated as reviewed labels. For segmented Stereo GEM, provide a `centroids_csv` with actual segmented centroids and an explicit coordinate scale in the measurement declaration. A transcript centroid is not a substitute for a segmentation centroid. Review the exported sheet and supply the completed file below.

```bash
vtx2 prepare path/to/square_008um context_v3.json prepared.json \
  --platform visium_hd --annotations existing_labels.csv --measurement measurement.json \
  --center 1200 1600 --max-cells 1000
```

The context must declare `schema_version: "3.0"`, `assay`/`platform: "visium_hd"`, `kernel: "lattice"` and spatial interactions. Example measurement:

```json
{
  "platform": "visium_hd",
  "unit": "binned_spot",
  "unit_size_um": 8,
  "microns_per_pixel": 0.25,
  "deconvolution_method": "replace with the reviewed method and reference",
  "purity_threshold": 0.7
}
```

The scale above is an illustrative declaration, not a platform default. CSV rows require `cell_id,cell_type,lineage`; bins additionally need a reviewed `purity` estimate. Units with absent/insufficient purity or undeclared deconvolution remain `unknown` and REST-only. They remain occupied; unknown is never equivalent to empty.

The lattice preserves bin IDs/XY. `occupancy=0` must be explicitly declared for an empty measured bin. A filtered matrix with only occupied bins offers no empty expansion sites by default. Lattice occupancy is an abstract quantity, not a cell count; ensembles reject mixing lattice and cell-agent units.
