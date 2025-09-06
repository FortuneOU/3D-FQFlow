
% ========================================================================
% Project: Ultrasound Flow Simulation & Reconstruction Framework
% ------------------------------------------------------------------------
% Description:
%   This project implements a complete pipeline for ultrasound flow imaging
%   simulation, beamforming-based reconstruction, image post-processing, 
%   and quantitative evaluation. It leverages the MUST (Matlab UltraSound Toolbox) 
%   for RF simulation with GPU acceleration and supports both 2-D and 3-D modes.
%
% Workflow Overview:
%   1) Define experiment parameters (probe, phantom, motion models).
%   2) Generate flow and tissue scatterers from CFD/phantom dataset.
%   3) Simulate ultrasound RF signals (FQSim2D/FQSim3D).
%   4) Perform Delay-and-Sum beamforming (FQBmFr2D/FQBmFr3D).
%   5) Apply SVD-based clutter filtering to separate flow/tissue signals.
%   6) Obtain Ultrasound Power Doppler Images (uPDI).
%   7) Quantitatively evaluate reconstruction accuracy via MSE, PSNR, SSIM, NCC.
%
% Features:
%   - Supports multiple phantoms (Renal, AiVsl, Gln).
%   - Supports different probes (L11-5v, Vermon).
%   - Motion models: static tissue, uniform flow, rotational field,
%     realistic optical-flow derived tissue motion.
%   - Highly scalable: GPU support, distributed computation, 
%     2D/3D volumetric reconstruction.
%   - Modular design with components: 
%        * alg/      : Core simulation & beamforming algorithms
%        * simcore/  : Scene, scatterers, reconstruction & processing classes
%        * utils/    : Helper utilities (I/O, transforms)
%
% Author: Qiang Fu
% Institution: Peking University
% Date: 2025-08-21
% ========================================================================


% Main script for simulation, reconstruction and post-processing
% -------------------------------------------------------------
% Adds paths, defines simulation/reconstruction parameters,
% runs simulation engine, beamforming reconstruction, 
% image processing, and final quantification.

% Add required toolboxes to path
addpath('E:\FQFlowNew\MUST');
addpath('.\simcore','.\alg','.\utils');

% ----- GLOBAL SIMULATION SETTINGS -----
globalParam.phantomCase = 'Renal';    % Phantom type: {Renal, AiVsl, Gln}
globalParam.probeCase   = 'Vermon';  % Probe type: {L11-5v, Vermon}
globalParam.MotionMode  = 0;         % 0: static tissue, 1: tissue moves, 2: vessel moves with tissue
globalParam.SimThreeDMode   = 1;     % 0: 2D simulation, 1: 3D simulation
globalParam.ReconThreeDMode = 1;     % 0: 2D reconstruction, 1: 3D
globalParam.velocityFieldMode = 0;   % Tissue motion field
                                     % 0: none, 1: uniform, 2: rotation, 3: real (from video)
globalParam.flowRCcoef  = 1/30;      % Reflection coefficient scaling factor

% Define reconstruction area
[xi,yi,zi] = meshgrid(linspace(-2e-2,2e-2,100), ...
                      linspace(0e-2,0e-2,100), ...
                      linspace(0e-2,4.5e-2,100));
globalParam.ReconArea.xi = xi;
globalParam.ReconArea.yi = yi;
globalParam.ReconArea.zi = zi;

% Define result folder
globalParam.resultFloder = '250830';

% ----- SIMULATION PHASE -----
sac = simcore.Scatterers(globalParam);     % Define scatterers
sac = sac.loadFlowScatters();              % Load flow scatterers
flowScatters = sac.flowScatters;

usParam = simcore.USParam(globalParam);    % Define probe parameters
scene   = simcore.Scene(globalParam);      % Define tissue scatterers & phantom
sim     = simcore.Simulator(flowScatters,scene,globalParam,usParam);

sim.runSimulation();                       % Run ultrasound simulation

% ----- RECONSTRUCTION PHASE -----
recon = simcore.Reconstructor(globalParam,usParam);
recon.runRecon();

% ----- IMAGE POST-PROCESSING -----
Imgp = simcore.ImageProcessor(globalParam);
Imgp.process2DImages();

% ----- QUANTIFICATION -----
qua = simcore.Quantification(globalParam);
qua.runQuantification();
