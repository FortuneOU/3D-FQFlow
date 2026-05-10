
# 🌊 FQ-FLOW: 3D-FQFlow Framework

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# 3D-FQFlow

Welcome to 3D-FQFlow (3D-Fully Quantitative Flow)! This is an open-source framework for quantitatively simulating blood flow and tissue motion for ultrafast power Doppler imaging (uPDI) . The tool is capable of simulating high-fidelity 3D vascular geometries, transient hemodynamics, and physiological tissue motion, while providing a highly optimized GPU-accelerated acoustic reconstruction module .

## 🧩 Module Overview

This framework consists of four interconnected core modules :

1. **Flow Simulation Module:** Employs computational fluid dynamics (CFD) to generate physiologically realistic blood flow kinematics. It integrates stochastic vascular generation, transient Navier-Stokes fluid dynamics (using the svFSI solver), and a custom high-precision Lagrangian particle tracking algorithm .
2. **Tissue Motion Simulation Module:** Generates perivascular tissue models and enables 3D motion simulation. Users can use parametric kinematic models or import real tissue motion fields extracted from clinical data using optical flow
3. **Ultrasound Simulation Module:** Implements a memory-efficient and GPU-accelerated frequency-domain ultrasound simulator (optimized PFIELD), supporting distributed parallel computing to rapidly generate large-scale 3D radio frequency (RF) imaging data containing millions of scatterers .
4. **Image Reconstruction and Processing Module:** Integrates an optimized 3D uPDI Delay-and-Sum (DAS) accelerated reconstruction algorithm, multiple image post-processing methods for B-mode/uPDI (including SVD filtering), and analysis tools for objective quantitative evaluation (MSE/PSNR/SSIM).

## 🚀 Basic Usage

Please follow these basic configuration steps to use the framework:

1. **Global Setup:** Defining a basic master configuration file (e.g., setting the spatial domain size and fundamental blood flow parameters).
2. **Environment Paths:** Specifying the local directory paths for the svFSI solver and other external tools within the master script.
3. **Particle Tracking Setup:** In the custom Python tracking module, configuring the model's flow inlet(s) and specifying the number of scatterers to be injected .
4. **Ultrasound & Motion Setup:** In the FQ-FLOW tissue motion simulator, defining the position of the imaging phantom, the ultrasound imaging parameters (e.g., transducer properties), and the macroscopic tissue motion parameters .

## 🎯 Minimal Working Example

To lower the barrier to entry for non-expert users, the framework provides an automated pipeline. A complete simulation can be executed with minimal manual intervention :

1. First, run `svvExample` to automatically generate and obtain a meshed fluid dynamics structural model.
2. Next, run `svFSI` to perform a complete CFD simulation on this meshed model (details for these two parts can be found in the official svv and svfsi documentation).
3. Then, call `particleTracker` to calculate and track the position coordinates of all scatterers (a set of pre-calculated coordinate examples is already provided in the code repository for direct use) .
4. Finally, run the `main_example.m` file in the `code` folder. The system will automatically read the coordinates and transducer parameters, and output the delay-and-sum (DAS) reconstructed multidimensional array, completing the entire calculation and imaging rendering process with one click .
---

## 🔗 External Dependencies
FQ-FLOW synergistically integrates several specialized solvers and libraries. To ensure full functionality, please install or configure the following dependencies based on the modules you intend to use:

### 1. Flow Simulation Module
*   **[svVascularize](https://github.com/SimVascular/svVascularize):** For stochastic constructive optimization and vascular network generation.
*   **[svFSI / svFSIplus](https://github.com/SimVascular/svFSI):** Core multiphysics finite element package for transient hemodynamic simulations.
*   **Python (3.x):** Required for automation scripts and managing the custom particle tracking module.
*   **[ParaView](https://www.paraview.org/):** Used for VTK data processing and extraction of flow field information.

### 2. Tissue Simulation Module
*   **MATLAB:** Used for generating and managing perivascular tissue motion.

### 3. Ultrasound Simulation Module
*   **MATLAB:** Core environment for acoustic field calculations.
*   **MATLAB Parallel Computing Toolbox:** **Required** for hardware acceleration to handle large-scale scatterer simulations efficiently.
*   **PFIELD (via [SIMUS](https://www.creatis.insa-lyon.fr/site7/en/SIMUS)):** Integrated component for simulating acoustic fields and transducer responses.

### 4. Image Processor
*   **MATLAB:** For Doppler signal reconstruction, filtering, and quantitative image analysis.

> [!NOTE]
> While the core Flow Simulation components (svFSI/svVascularize) are open-source, the Ultrasound and Tissue modules require a **MATLAB license**. Please ensure you have the necessary subscriptions, including the *Parallel Computing Toolbox*, for full performance.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
