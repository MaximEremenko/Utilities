# Diffuse Data Format Converter

A browser-based converter for 3-D single-crystal diffuse scattering data.
It translates the same intensity grid between the file formats used by
RMCProfile, DISCUS, Yell, Meerkat, and Scatty, entirely inside the browser.
No installation, no server, and no upload: files never leave your computer.

## Start Here

- [Live Tool](Format_Converter.html)
- Or open `Format_Converter.html` locally in any modern browser
  (double-click works; no web server is needed).

## Supported Formats

| Format | Read | Write | Typical producer |
| --- | --- | --- | --- |
| Unified data format (HDF5) | yes | yes | RMCProfile, DiffuseCode/DISCUS |
| Yell 1.0 (HDF5) | yes | yes | DISCUS, Yell, Meerkat |
| RMCProfile old text format (.dat) | yes | yes | RMCProfile Diffuse3D |
| Scatty VTK (.vtk) | yes | yes | Scatty; also loads in ParaView |

Any format can be converted to any other. The input format is detected
automatically from the file content.

### Format notes

- **Unified data format**: both the `/scattering/data` compatibility layout
  and the RMCProfile/DISCUS `/entry/data` layout are read; written files
  contain both. New output follows the current common contract: the mandatory
  crystal-metadata audit datasets are written under `/entry/data`, with
  `audit_conform_dict_name = Disorder unified data`,
  `data_type_axes = hkl`, and `data_type_number = real`. Readers accept
  the older `Disorder scattering` name and infer missing axes metadata as
  RMCProfile does for legacy files, but reject a structure dictionary or an
  unsupported axes/number type instead of silently treating it as HKL data.
  Three-dimensional `Q` axes are converted to HKL using the real parent
  cell; scalar or direct-space axis types are outside this diffuse converter.
  Optional crystal flags do not control diffuse-grid interpretation; the
  standard cell and symmetry fields are always written.
- **Yell 1.0**: the flat layout with `/data`, `/lower_limits`, `/step_sizes`,
  and `/unit_cell` datasets. Files that store the unit metric
  (`unit_cell = 1 1 1 90 90 90`) are interpreted as hkl grids in reciprocal
  lattice units of the parent cell; written files always store the real cell
  when one is known.
- **RMCProfile old text format**: an `npoints nsec` header followed by
  `i j k qx qy qz intensity` rows. Q is cartesian in 1/Angstrom with the
  2*pi convention, `q = 2*pi * B * hkl`, where `B` is the reciprocal basis
  of the parent cell. Files with several symmetry sections are read using
  the section 1 coordinates.
- **Scatty VTK**: ASCII `STRUCTURED_POINTS` as written by Scatty. `ORIGIN`
  and `SPACING` are cartesian Q in 1/Angstrom (2*pi convention, the same as
  the old text format); values run x fastest, z slowest. Writing VTK
  requires the Q grid to be axis-aligned and ascending, because
  `STRUCTURED_POINTS` cannot express rotated or non-orthogonal grids
  (Scatty itself has the same restriction).

## Usage

1. **Data file**: select the diffuse data file (`.h5`, `.dat`, or `.vtk`).
   The format, grid size, and stored cell are reported in the log.
2. **Unit cell**: some conversions need the parent (crystallographic) unit
   cell, because the text and VTK formats do not store one, and some Yell
   files store only the unit metric. Provide either
   - a structure file: `.rmc6f` or unified structure `.h5`. For `.rmc6f`,
     the parent cell is the supercell divided by its dimensions. In the
     unified structure contract, `unit_cell_lengths` already stores the
     basic/parent cell and is used directly; or
   - the parent cell typed in directly (a, b, c in Angstrom; alpha, beta,
     gamma in degrees).
   When the data file already stores a real cell (for example RMCProfile
   `_calc.h5` output), it is used automatically and this section can be left
   empty.
3. **Output**: pick the target format, optionally set the radiation metadata
   for HDF5 output, and press "Convert and download".

### When is the unit cell required?

| Input | Cell needed? |
| --- | --- |
| Unified `.h5` with a real cell | no |
| Yell `.h5` with a real cell | no |
| Yell `.h5` with unit metric | yes (for `.dat`/`.vtk` output; passed through for HDF5 output) |
| Old text `.dat` | yes, always |
| Scatty `.vtk` | yes, always |

## Examples

The `Examples/` folder contains one small synthetic dataset (cubic parent
cell a = 5.63 Angstrom, 5 x 5 x 5 grid, hkl from -1 to 1 in steps of 0.5)
written in every supported format, plus a matching structure file:

- `example_structure.rmc6f` - structure with a 2 x 2 x 2 supercell
- `example_unified.h5` - unified data format
- `example_yell.h5` - Yell 1.0
- `example_diffuse3d.dat` - RMCProfile old text format
- `example_scatty.vtk` - Scatty VTK

Load any of the data files, add `example_structure.rmc6f` where a cell is
required, and convert in any direction; the intensity values are identical
in all files, so results are easy to compare.

## Included Surface

- `Format_Converter.html`
  - Main browser app.
- `js/converter.js`
  - Format readers/writers and the cell/reciprocal-space math. Plain
    JavaScript, also loadable from Node.js for testing.
- `js/h5wasm.js`
  - Vendored [h5wasm](https://github.com/usnistgov/h5wasm) 0.10.3 bundle
    (the HDF5 library compiled to WebAssembly, self-contained). License in
    `js/h5wasm-LICENSE.txt`.
- `Examples/`
  - The example dataset described above.

## Validation

The conversion core was verified against real files from each producer:

- Yell/DISCUS: a DISCUS-written Yell 1.0 file converts to the RMCProfile old
  text format reproducing an independently produced reference `.dat` for the
  same dataset at the reference file's own precision.
- Scatty: parsing a Scatty `.vtk` reproduces the intensities of Scatty's
  paired `*_list.txt` hkl output exactly, and the recovered hkl grid matches
  it to the 6-decimal precision of the VTK header.
- RMCProfile: HDF5 files written by the converter follow the same dataset
  layouts, mandatory audit/type metadata, axis ordering, and fixed-length
  string types as RMCProfile's Fortran implementation
  (`unified_config/unified_hdf5_io.f90`), and
  round-trip conversions through every format are exact at double precision.

## Command-Line Equivalents

RMCProfile (built with `RMC_ENABLE_HDF5=ON`) provides Fortran tools for the
unified/old-text conversions: `unified_to_diffuse3d` and
`diffuse3d_to_unified`.
