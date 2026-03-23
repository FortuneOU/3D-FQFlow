% main.m for Sequence Image processing

% Add required toolboxes to path
addpath('..\MUST');
addpath('.\simcore','.\alg','.\utils');

% ----- GLOBAL SIMULATION SETTINGS -----
globalParam.phantomCase = 'F20Trace';    % Phantom type: {Renal, AiVsl, Gln}
globalParam.probeCase   = 'Vermon';  % Probe type: {L11-5v, Vermon，Vermon-3}
globalParam.MotionMode  = 0;         % 0: static tissue, 1: tissue moves, 2: vessel moves with tissue %10 no tissue
globalParam.SimThreeDMode   = 1;     % 0: 2D simulation, 1: 3D simulation
globalParam.ReconThreeDMode = 1;     % 0: 2D reconstruction, 1: 3D
globalParam.velocityFieldMode = 0;   % Tissue motion field
                                     % 0: none, 1: uniform, 2: rotation, 3: real (from video)
globalParam.flowRCcoef  = 0.1;      % Reflection coefficient scaling factor

% Define reconstruction area
% [xi,yi,zi] = meshgrid(linspace(-2e-2,2e-2,200), ...
%                       linspace(0,0,1), ...
%                       linspace(0.5e-2,5.5e-2,250));
[xi,yi,zi] = meshgrid(linspace(-0.5e-2,0.5e-2,100), ...
                      linspace(-0.5e-2,0.5e-2,100), ...
                      linspace(0.5e-2,1.5e-2,100));
globalParam.ReconArea.xi = xi;
globalParam.ReconArea.yi = yi;
globalParam.ReconArea.zi = zi;


% Define result folder
globalParam.resultFloder = '260317_F20_linear';


% ----- SIMULATION PHASE -----
 usParam = simcore.USParam(globalParam);    % Define probe parameters
scene   = simcore.Scene(globalParam);      % Define tissue scatterers & phantom
% 
sac = simcore.Scatterers(globalParam);     % Define scatterers
sac = sac.loadFlowScattersFromMat4Seq();
flowScatters = sac.flowScatters;

sim = simcore.Simulator(flowScatters,scene,globalParam,usParam);
sim.runSimulation();                       % Run ultrasound simulation

% ----- RECONSTRUCTION PHASE -----
recon = simcore.Reconstructor(globalParam,usParam);
recon.runRecon();


% % ----- IMAGE POST-PROCESSING -----
Imgp = simcore.ImageProcessor(globalParam,usParam);
%Imgp.process2DImages4BanduPDI();
%Imgp.process2DImages4C();
%Imgp.process2DImages4PW();

Imgp.process3DImages();
