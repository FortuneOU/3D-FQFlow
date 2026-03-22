
# 🌊 FQ-FLOW: 3D-FQFlow Framework

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**FQ-FLOW (3D-FQFlow)** is a comprehensive framework that synergistically integrates multiple open-source hemodynamic simulation tools with highly optimized ultrasound simulation and imaging reconstruction algorithms. It provides a full-pipeline solution from stochastic vascular generation and fluid dynamics to realistic 3D ultrasound radio-frequency (RF) signal synthesis and image reconstruction.

---

## 🏗 Architecture & Modules

The framework comprises four interconnected simulation modules to model blood vessels, perivascular tissues, RF signal simulation, and image reconstruction:

### 1. Flow Simulation Module
Generates physiologically realistic blood flow kinematics to serve as the ground truth for ultrasound scattering modeling.
*   **Stochastic Vascular Generation:** Utilizes constrained constructive optimization to synthesize plausible microvascular architectures based on Murray's law.
*   **Hemodynamic Simulation:** Solves transient Navier-Stokes equations on unstructured meshes using a multiphysics finite element solver to capture complex, pulsatile flow profiles (e.g., systolic/diastolic variations).
*   **Lagrangian Tracking:** Features a custom Python-based discrete-phase Lagrangian tracking algorithm using 2nd-order Runge-Kutta integration to accurately track ultrasound scatterers (RBCs/microbubbles). Includes a dynamic resetting mechanism for stagnant particles to ensure stable scatterer concentration.

### 2. Tissue Motion Simulation Module
Simulates perivascular tissue models and realistic 3D motion, which is crucial for evaluating clutter filtering in ultrafast Doppler imaging (uPDI).
*   **Tissue Microstructure:** Employs randomly distributed point scatterers (~83 scatterers/mm³) based on a Rayleigh distribution.
*   **Kinematic Integration:** Accurately couples perivascular tissue motion with intravascular hemodynamics via rigid-body kinematic superposition.
*   **Real Data Driven:** Capable of extracting motion fields directly from real clinical ultrasound images using optical flow methods.

### 3. Ultrasound Simulation Module
A highly optimized, parallelized computing architecture for generating large-scale 3D RF data based on the weak scattering assumption.
*   **Static Tissue Optimization:** Separates blood flow scatterers from static tissue scatterers. Tissue RF echoes are computed only once and reused, dramatically reducing redundant computations.
*   **GPU-Accelerated PFIELD:** Offloads intensive frequency-domain matrix multiplications to the GPU and uses a dynamic distributed memory strategy. It successfully simulates $1 \times 10^6$ scatterers in ~4,117 seconds with only 1.9 GB of peak memory consumption, avoiding the Out-Of-Memory (OOM) issues of existing simulators.

### 4. Image Reconstruction and Processing Module
Integrates optimized algorithms and evaluation tools for post-processing.
*   **Accelerated uPDI Reconstruction:** Uses distributed Delay-and-Sum (DAS) processing with precomputed matrices, significantly slashing multi-frame reconstruction times (e.g., 30 frames reconstructed in 52 mins vs. 981 mins in conventional methods).
*   **Post-Processing:** Supports both B-mode and uPDI imaging (2D/3D) with SVD filtering and logarithmic compression.
*   **Quantitative Analysis:** Built-in evaluation tools using Mean Squared Error (MSE), Peak Signal-to-Noise Ratio (PSNR), and Structural Similarity Index (SSIM).

---

## 🚀 Minimal Working Example (MWE) & Automation

FQ-FLOW is designed to act as an automated pipeline, lowering the technical barrier for non-expert users. You act as the *experimental designer* through high-level configurations:

1.  **Configuration:** Define a master configuration file specifying the spatial domain, flow parameters, paths to solvers, flow inlets, scatterer populations, and imaging/transducer properties.
2.  **Automated Execution:** 
    *   The pipeline autonomously invokes the vascular generator to build the stochastic network and 3D mesh.
    *   The mesh is passed to the hemodynamic solver to compute multi-time-step 3D flow fields.
    *   Python tracking programs resolve Lagrangian trajectories across all time steps.
    *   Scatterer positions and tissue motion are fed into the acoustic simulator to render the ultrasound sequence.

By simply modifying the configuration parameters, complex 3D simulations can be reproduced without deep expertise in manual software coupling.

---

## 💻 How to Use

### Prerequisites
Ensure you have Python 3.8+ installed along with the required libraries:
```bash
pip install numpy scipy pyvista pyyaml
```

### Quick Start
1.  **Clone the repository:**
    ```bash
    git clone https://github.com/yourusername/FQ-FLOW.git
    cd FQ-FLOW
    ```
2.  **Configure the Master Settings:**
    Open `config/master_config.yaml` and set your spatial domain, paths to external solvers, transducer properties, and simulation parameters.
3.  **Run the Pipeline:**
    Execute the automated pipeline script:
    ```bash
    python run_pipeline.py --config config/master_config.yaml
    ```
4.  **View Results:**
    RF signals and reconstructed images (B-mode/uPDI) will be saved in the `results/` directory. You can use the provided evaluation scripts to run MSE/PSNR/SSIM analysis.

---

## 🔗 External Dependencies

FQ-FLOW synergistically integrates world-class open-source solvers. Please ensure you have them installed or compiled according to your system requirements:

*   **[svVascularize](https://github.com/SimVascular/svVascularize):** Used in the Flow Simulation Module for stochastic constructive optimization and vascular network generation.
*   **[svFSI / svFSIplus](https://github.com/SimVascular/svFSI):** A highly capable multiphysics finite element package used for transient hemodynamic simulations.

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
