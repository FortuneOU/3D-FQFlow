# 3D-FQFlow
An Open-Source Platform for Fully Quantitative Ultrasound Flow and Tissue Motion Simulation

Flow Simulator
This project provides a comprehensive flow simulation pipeline for modeling, reconstructing, and visualizing three-dimensional vascular structures and associated hemodynamic fields. It leverages cutting-edge computational and visualization tools to generate physiologically realistic vascular networks, simulate blood flow, analyze particle dynamics, and produce publication-quality visualizations. The pipeline features the following core components:

1. Deep Learning-based 3D Vascular Structure Generator
At the heart of the simulator lies a stochastic parametric L-system-based generator for synthetic vascular structures. This module simulates fractal vascular growth through iterative probabilistic rewriting rules, constrained by biomechanical and anatomical parameters:

Fractal Growth: Models vascular bifurcation, segment length, and diameter variation according to biological rules (e.g., Murray's law, realistic bifurcation angles).
Stochasticity: Introduces natural anatomical variability and enables the simulation of vascular anomalies (e.g., stenosis, aneurysms).
3D Constraints: Accepts user-defined boundary surfaces to constrain vessel growth within specified anatomical domains.
Realistic Imaging: Generates simulated CT/MRI intensity profiles mimicking angiographic characteristics.

Flow/Structuer Generator/generate_vessel_network.py

2. SimVascular Integration
The generated vascular skeletons are readily imported into SimVascular , a widely used open-source platform for vascular modeling and simulation:

Centerline and Cross-section Extraction: Import synthetic structures to reconstruct full 3D vascular models.
Mesh Generation: Use TetGen to build unstructured tetrahedral meshes supporting complex geometries.
CFD Simulation: Employ svSolver (Finite Element Method) to solve the unsteady Navier–Stokes equations for blood flow.
https://simvascular.github.io/

3. ParaView Visualization
Results from SimVascular are visualized and post-processed in ParaView, a powerful open-source scientific visualization platform:

Data Processing: Handle structured/unstructured meshes and volumetric data with modular filters.
Field Visualization: Render velocity and pressure fields, extract slices, and create streamlines for in-depth flow analysis.
Parallelism: Large datasets can be processed efficiently thanks to ParaView's distributed architecture.
https://www.paraview.org/


4. PROTEUS Particle Trajectory Analysis
PROTEUS extends the flow simulator with a streamline-based model of microbubble or particle transport:

Streamline Integration: Computed from pre-solved CFD velocity fields using adaptive step-size control to ensure numerical accuracy.
Boundary Handling: Particles exiting the domain are probabilistically re-initialized at inflow boundaries according to velocity-weighted sampling.
Physically Realistic: Simulates complex 3D transport dynamics of scatterers under physiological flow, providing a platform for contrast agent studies.
https://github.com/PROTEUS-SIM/PROTEUS


Example Workflow
Generate: Use the 3D Vascular Structure Generator to create a synthetic vascular network.
Simulate: Import the network into SimVascular, generate a mesh, and run CFD simulations.
Visualize: Open simulation results in ParaView for field analysis and visualization.
Analyze: Apply PROTEUS to study the behavior of intravascular particles along flow streamlines.
