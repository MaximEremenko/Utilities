# Utilities Repository

Welcome to the Utilities repository by Maksim Eremenko. This repository hosts
browser-based tools for working with RMC6f configurations, diffuse scattering,
and related structural-analysis workflows.

## Start Here

Public-facing documentation for the diffuse-scattering utility now lives in the
GitHub Pages docs surface:

- [Diffuse Scattering Guide](https://maximeremenko.github.io/Utilities/docs/diffuse_scattering.html)
- [Diffuse Scattering Examples](https://maximeremenko.github.io/Utilities/docs/diffuse_scattering_examples.html)
- [Diffuse Scattering Theory](https://maximeremenko.github.io/Utilities/docs/diffuse_scattering_theory.html)
- [Diffuse Scattering Supplementary](https://maximeremenko.github.io/Utilities/docs/diffuse_scattering_supplementary.html)
- [Diffuse Scattering Bibliography](https://maximeremenko.github.io/Utilities/docs/diffuse_scattering_bibliography.html)
- [Diffuse Scattering Troubleshooting](https://maximeremenko.github.io/Utilities/docs/diffuse_scattering_troubleshooting.html)
- [Utilities Docs Home](https://maximeremenko.github.io/Utilities/docs/index.html)

## Live Demos

- **Utilities Home:** [Utilities Home](https://maximeremenko.github.io/Utilities/)
- **Diffuse Scattering Tool:** [RMC6f Diffuse Scattering Tool](https://maximeremenko.github.io/Utilities/RMCProfileUtilities/Diffuse_Scattering/Diffuse_Scattering_RMC6f.html)
- **Select Subvolume App:** [RMC6f Select Subvolume App](https://maximeremenko.github.io/Utilities/RMCProfileUtilities/Select_Subvolume/RMC6f_Select_Subvolume_App.html)
- **PCA KDE Tool:** [PCA KDE Tool](https://maximeremenko.github.io/Utilities/RMCProfileUtilities/PCA_KDE/PCA_KDE_rmcdisplacements.html)
- **PCA SDE Tool:** [PCA SDE Tool](https://maximeremenko.github.io/Utilities/RMCProfileUtilities/PCA_SDE/PCA_SDE_rmcdisplacements.html)
- **Slice Diffuse 3D Tool:** [RMC Diffuse 3DSlice Tool](https://maximeremenko.github.io/Utilities/RMCProfileUtilities/RMCDiffuse3DSlice/SliceDiffuse3DRMC.html)
- **Background Remover:** [Background Remover Tool](https://maximeremenko.github.io/Utilities/RMCProfileUtilities/Background_Remover/Background_Remover.html)
- **Format Converter:** [Diffuse Data Format Converter](https://maximeremenko.github.io/Utilities/RMCProfileUtilities/Format_Converter/Format_Converter.html)

## Repository Structure

- `docs/`
  - GitHub Pages documentation entry point.
- `RMCProfileUtilities/`
  - Source tree for the browser tools.
- `RMCProfileUtilities/Diffuse_Scattering/`
  - In-browser diffuse-scattering calculator for `.rmc6f` files.
- `RMCProfileUtilities/Diffuse_Scattering/js/`
  - Reader, scattering, and diffuse-calculation modules.
- `RMCProfileUtilities/Format_Converter/`
  - In-browser converter for 3-D diffuse scattering data formats.

## Diffuse Scattering Tool

The diffuse-scattering tool computes and visualizes `I(h,k,l)` directly in the
browser from `.rmc6f` input files. It provides:

- 3D slice-plane visualization
- 2D slice-map visualization
- 3D isosurface rendering
- axis slices, normal-plane slices, and volume-average slices
- neutron, X-ray, and electron scattering options
- CPU and WebGPU computation paths
- export to `.dat`, `.json`, `.vtk`, `.cube`, slice CSV/SVG, and Plotly HTML/SVG

For the detailed user guide, examples, and theory pages, use the docs links at
the top of this README.

## Format Converter

The format converter translates 3-D diffuse scattering data between the file
formats used by RMCProfile, DISCUS, Yell, Meerkat, and Scatty, entirely in the
browser:

- unified data format HDF5 (`/scattering/data` and `/entry/data` layouts)
- Yell 1.0 HDF5
- RMCProfile old text format (`.dat`)
- Scatty VTK (`.vtk`)

Conversions that need the parent unit cell accept an `.rmc6f` or unified
structure file, or a manually entered cell. See
[Format_Converter/README.md](RMCProfileUtilities/Format_Converter/README.md)
for details and shipped examples.

## Supplementary Material

Backend provenance, scattering-table families, and reproducibility notes are collected in:

- [Diffuse Scattering Supplementary](https://maximeremenko.github.io/Utilities/docs/diffuse_scattering_supplementary.html)

## Contribution

Feel free to open issues or submit pull requests for improvements and bug fixes.

## License

This repository is distributed under Apache-2.0. See [LICENSE](LICENSE).
