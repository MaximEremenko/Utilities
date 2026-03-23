# RMC6f Diffuse Scattering Tool

This directory contains the browser-based diffuse-scattering calculator for
`.rmc6f` configurations.

The tool computes diffuse intensity `I(h,k,l)` in reciprocal space, visualizes
3D and 2D slices, renders 3D isosurfaces, and exports the resulting data.

## Start Here

- [User Guide](../../docs/diffuse_scattering.html)
- [Examples](../../docs/diffuse_scattering_examples.html)
- [Theory](../../docs/diffuse_scattering_theory.html)
- [Troubleshooting](../../docs/diffuse_scattering_troubleshooting.html)
- [Live Tool](Diffuse_Scattering_RMC6f.html)

## Included Surface

- `Diffuse_Scattering_RMC6f.html`
  - Main browser app.
- `js/`
  - Modular source files for RMC6f parsing, scattering coefficients, and
    diffuse-calculation logic.
- `Examples/`
  - Shipped `.rmc6f` profiles for external users of the Utilities docs surface.

## Example Profiles Highlighted In The Docs

- `LiFeO2.rmc6f`
  - Chemical-order benchmark with diffuse manifold plus `1/2(111)` condensation.
- `CaTiO3.rmc6f`
  - Displacement benchmark with overlapping rod-like and breathing-related
    diffuse features.
- `PMN_300k.rmc6f`
  - Relaxor benchmark for anisotropic diffuse features and complex slice
    exploration.

## Scope

This tool is a forward diffuse-scattering calculator and visualization surface.
It does not perform the MOSAIC inverse reconstruction workflow; the MOSAIC paper
is referenced in the docs for scientific benchmark context only.
