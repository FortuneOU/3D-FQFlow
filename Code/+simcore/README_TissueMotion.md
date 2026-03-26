## 🫀 Tissue Motion Simulation Module

In practical clinical applications of ultrafast Power Doppler Imaging (uPDI), perivascular tissue motion induced by respiration and cardiac pulsation generates clutter signals (i.e., "flash noise") that are typically 20-30 dB stronger than the blood flow signals. Therefore, accurately simulating tissue motion is crucial for developing and validating advanced clinical clutter filtering algorithms.

### 📐 Basic Principles & Core Formulas

1. **Microstructure Characterization:** Tissue microstructure is modeled as randomly distributed point scatterers (with a density of approximately 83 scatterers/mm³). Their reflection coefficients are sampled from a Rayleigh distribution with a mean value of 1. The overall acoustic field response is then calculated using the principle of linear superposition.
2. **Kinematic Integration of Tissue and Blood Flow:** The module determines the absolute spatial trajectories of intravascular scatterers through rigorous kinematic superposition. The core physical formula can be expressed as:
   `V_total(x,y,z,t) = V_hemo(x,y,z,t) + V_tissue(x,y,z,t)`
   This means that the instantaneous absolute velocity of a blood scatterer (`V_total`) is the vector sum of the local Eulerian hemodynamic velocity driven by the cardiac cycle (`V_hemo`) and the global tissue motion vector at that specific coordinate (`V_tissue`). This rigid-body kinematic integration ensures physical continuity between the internal blood flow and the moving vascular boundaries.
3. **Realistic Motion Field Extraction:** In addition to defining parametric motions, this module supports the use of the Optical Flow method to directly extract complex, periodic tissue motion fields from real clinical ultrasound image data.

### 💻 Core Code Parsing (`Simulator` Class)

The configured `Simulator` class acts as the core engine driving the interaction between the physical processes described above and the underlying ultrasound algorithms. Its core workflow is divided into two phases: Trajectory Calculation and RF Synthesis:

1. **Motion Field Generation & Trajectory Integration (`velocityFieldMode`):**
   In the constructor, the system extracts the time step based on the ultrasound pulse repetition frequency (`dt = 1/usParam.PRF`). By integrating `points = points + velocities * dt`, the system calculates the tissue trajectory for each frame and saves it locally. The code natively supports three types of motion fields:
   - **Mode 1 (Uniform Field):** Sets a fixed spatial translation velocity (e.g., `Vx = 2mm/s, Vz = 4mm/s`).
   - **Mode 2 (Rotational Field):** Constructs a rotational velocity matrix around a specific center point.
   - **Mode 3 (Realistic Field):** Utilizes the `opticalFlowHS` algorithm to process imported clinical ultrasound videos (like the kidney data in the example). It then uses a 2D median filter (`medfilt2`) to smooth the extracted motion field and eliminate noise spikes.

2. **RF Synthesis & Computational Optimization (`runSimulation`):**
   - **Motion Coupling (`MotionMode`):** The engine determines the merging strategy based on user settings. For tissue-only motion, the moved `TissuePoints` are simply concatenated with the current blood flow frame. If the vessels move together with the tissue, the code applies a frame-cumulative displacement offset (`deltaD`) to the blood scatterers.
   - **Probe Field-of-View Cropping:** For different ultrasound probes (e.g., L11-5v or Vermon matrix arrays), the code automatically crops and removes scatterers that fall outside the physical detection region `(x, y, z)` before entering the acoustic computation phase. This significantly reduces GPU memory pressure.
   - **Ultra-Fast Rendering Strategy for Static Tissue:** As shown in the `needSimPhantom` logic, when there is no tissue motion (static condition), the module separates the RF echoes of the stationary extravascular tissue from the blood flow. The tissue echo is computed only once in the very first frame (`PhantomRF`). For all subsequent frames, the system only computes the time-varying pure blood flow RF, pads the matrices to match dimensions, and directly adds them together. This optimization strategy eliminates massive amounts of redundant inter-frame computations and is the key to enabling highly efficient large-scale 3D simulations.
