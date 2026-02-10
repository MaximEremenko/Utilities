# Utilities Repository

Welcome to the Utilities Repository by Maxim Eremenko. This repository contains a collection of useful utilities and analysis tools for working with RMC6f files.

## Live Demos

Access live demos via GitHub Pages:

- **Home Page:** [Utilities Home](https://maximeremenko.github.io/Utilities/)
- **Select Subvolume App:** [RMC6f Select Subvolume App](https://maximeremenko.github.io/Utilities/RMCProfileUtilities/Select_Subvolume/RMC6f_Select_Subvolume_App.html)
- **PCA KDE Tool:** [PCA KDE Tool](https://maximeremenko.github.io/Utilities/RMCProfileUtilities/PCA_KDE/PCA_KDE_rmcdisplacements.html)  
  *(Tool for performing Kernel Density Estimation as part of the RMC analysis)*
- **PCA SDE Tool:** [PCA SDE Tool](https://maximeremenko.github.io/Utilities/RMCProfileUtilities/PCA_SDE/PCA_SDE_rmcdisplacements.html)  
  *(Tool for performing Surface Density Estimation as part of the RMC analysis)*
- **Diffuse Scattering Tool (Client-side):** [RMC6f Diffuse Scattering Tool](https://maximeremenko.github.io/Utilities/RMCProfileUtilities/Diffuse_Scattering/Diffuse_Scattering_RMC6f.html)  
  *(Compute/visualize/save diffuse scattering I(h,k,l) in-browser from RMC6f input using Type-3 NUFFT CPU/WebGPU backends and selectable neutron/X-ray/electron scattering models)*
- **Slice Diffuse 3D RMCProfile Tool:** [RMC Diffuse 3DSlice Tool](https://maximeremenko.github.io/Utilities/RMCProfileUtilities/RMCDiffuse3DSlice/SliceDiffuse3DRMC.html)
- **Background Remover:** [Background Remover Tool](https://maximeremenko.github.io/Utilities/RMCProfileUtilities/Background_Remover/Background_Remover.html)  
  *(Interactive tool to remove white or near-white backgrounds from PNG, JPEG, and TIFF images. Useful for preprocessing images in RMC analysis and visualization workflows.)*

## Repository Structure

- **Utilities/**  
  The root of the repository.

- **Utilities/docs/**  
  GitHub Pages documentation entry point.

- **Utilities/RMCProfileUtilities/**  
  Contains tools related to RMC profile analysis.
  - **Select_Subvolume/** - Source code for the RMC6f Select Subvolume App.
  - **PCA_KDE/** - PCA KDE analysis tools for RMC6f files.
  - **PCA_SDE/** - PCA SDE analysis tools for RMC6f files.
  - **Diffuse_Scattering/** - In-browser diffuse scattering calculator for RMC6f files.
  - **Diffuse_Scattering/js/** - Modular scripts: RMC6f reader, diffuse core math, amplitude engine, and scattering coefficient sources/database.
  - **RMCDiffuse3DSlice/** - Utility to slice RMCProfile diffuse 3D DAT files.
  - **Background_Remover/** - Utility to remove near-white image backgrounds.

## Analysis Tools for RMC6f Files

This repository includes several utilities designed to help analyze RMC6f files:

- **Select Subvolume App:**  
  A tool for selecting subvolumes within RMC6f files.

- **PCA KDE Tool:**  
  Utilizes Kernel Density Estimation (KDE) as part of the RMC analysis.

- **PCA SDE Tool:**  
  Utilizes Surface Density Estimation (SDE) as part of the RMC analysis.

- **Diffuse Scattering Tool (Client-side):**  
  Computes diffuse scattering from RMC6f configurations in the browser, visualizes 3D + 2D slices, and saves `.dat`/`.json` outputs. Supports neutron/X-ray/electron scattering options and CPU/WebGPU NUFFT paths.

- **Slice Diffuse 3D RMCProfile Tool:**  
  Tool to make slices of RMCProfile Diffuse 3D dat files.

- **Background Remover:**  
  Interactive app for removing white or near-white backgrounds from images. Supports PNG, JPEG, and TIFF formats. Useful for preparing images for presentations or further analysis.

## How to Contribute

Feel free to open issues or submit pull requests if you have improvements or bug fixes.

---

For any questions or support, please contact Maxim Eremenko.
