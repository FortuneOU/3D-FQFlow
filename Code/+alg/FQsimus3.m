function [RF,param,RFspectrum] = FQsimus3(varargin)
%FQsimus3 - Simulation of ultrasound RF signals for a planar 2-D array
%   [RF,PARAM] = FQsimus3(X,Y,Z,RC,DELAYS,PARAM) simulates the RF
%   radio-frequency signals generated when a planar 2D array transmits
%   specified wavefronts through a scattering medium.
%
%   INPUTS:
%     - X,Y,Z : scatterer coordinates [m]
%     - RC    : reflection coefficients (same size as X)
%     - DELAYS: transmit delay laws [s]
%     - PARAM : transducer parameters structure
%
%   OUTPUTS:
%     - RF        : simulated RF data matrix (#samples × #channels)
%     - param     : updated parameter structure
%     - RFspectrum: RF frequency-domain signals
%
%   Notes:
%     * Supports multi-plane transmission (MPT) via DELAYS matrix
%     * Uses PFIELD3 during TX/RX modeling
%     * Implements GPU acceleration & distributed computing (multi-chunk)
%     * Sampling frequency param.fs defaults to 4 × fc (if not set)
%
%   Reference:
%     SIMUS3 is part of the MUST (Matlab UltraSound Toolbox).
% -----------------------------------------------------------------------
% Author: Damien Garcia (MUST), 2022. Modified version with GPU support.
% -----------------------------------------------------------------------

if nargin==0
    if nargout>0
        [RF,param] = RunTheExample;
    else
        RunTheExample;
    end
    return
end
narginchk(6,7)
nargoutchk(0,3)

% ------------ INPUT PARSING ----------------
x = varargin{1};
switch nargin
    case 6
        y = varargin{2}; z = varargin{3}; RC = varargin{4};
        delaysTX = varargin{5}; param = varargin{6};
        options = [];
    otherwise
        y = varargin{2}; z = varargin{3}; RC = varargin{4};
        delaysTX = varargin{5}; param = varargin{6};
        options = varargin{7};
end

% Validate scatterer arrays
assert(isequal(size(x),size(y),size(z),size(RC)),...
    'X, Y, Z, and RC must be identical size')
if isempty(x)
    RF = []; RFspectrum = []; return
end

% Normalize field names
param   = IgnoreCaseInFieldNames(param);
options = IgnoreCaseInFieldNames(options);
options.CallFun = 'simus3';

% Default waitbar
if ~isfield(options,'WaitBar'); options.WaitBar = true; end
assert(islogical(options.WaitBar),'OPTIONS.WaitBar must be logical')

% Default parallel setting
if ~isfield(options,'ParPool'); options.ParPool = false; end

% Sampling frequency
if ~isfield(param,'fs')
    param.fs = 4*param.fc; 
end
assert(param.fs>=4*param.fc,'PARAM.fs must be >= 4*fc')

NumberOfElements = size(param.elements,2);

% Reception delays
if ~isfield(param,'RXdelay')
    param.RXdelay = zeros(1,NumberOfElements);
else
    assert(numel(param.RXdelay)==NumberOfElements,...
        'PARAM.RXdelay length mismatch with #elements')
end

% Spectrum threshold
if ~isfield(options,'dBThresh')
    options.dBThresh = -100;
end

% Frequency step factor
if ~isfield(options,'FrequencyStep')
    options.FrequencyStep = 1;
end
if options.FrequencyStep>1
    warning('OPTIONS.FrequencyStep > 1: aliasing risk')
end

% ------------ COMPUTATION PREPARATION ----------------
xe = param.elements(1,:);
ye = param.elements(2,:);
d2 = (x(:)-xe).^2 + (y(:)-ye).^2 + z(:).^2;
maxD = max(sqrt(d2(:))); 
clear d2

% Add pulse length to max depth
[~,tp] = getpulse(param,2);
maxD = maxD + tp(end)*param.c;

% Parallel pool (if enabled)
options.ParPool = options.ParPool & license('test','Distrib_Computing_Toolbox');
if options.ParPool
    pool = gcp;
    NW = pool.NumWorkers; Nx = numel(x);
    dim1Dist = [ones(1,NW-1)*floor(Nx/NW) Nx-(NW-1)*floor(Nx/NW)];
    x  = mat2cell(x(:),dim1Dist,1);
    y  = mat2cell(y(:),dim1Dist,1);
    z  = mat2cell(z(:),dim1Dist,1);
    RC = mat2cell(RC(:),dim1Dist,1);
end

% Frequency sampling to avoid aliasing in time domain
df = 1/2/(2*maxD/param.c + max(delaysTX+param.RXdelay,[],'all'));
df = df*options.FrequencyStep;
Nf = 2*ceil(param.fc/df)+1;
RFspectrum = zeros(Nf,NumberOfElements);
options.FrequencyStep = df;

% ------------ CALL PFIELD3 TO COMPUTE RF SPECTRA ------------
if options.ParPool
    options.WaitBar = false;
    spmd(NW)
        options.RC = RC{labindex};
        [~,~,RFsp,idx] = alg.FQpfield3(x{labindex},y{labindex},z{labindex},...
            delaysTX,param,options);
    end
    for k = 1:NW
        RFspectrum(idx{k},:) = RFspectrum(idx{k},:) + RFsp{k};
    end
else
    % Without pool -> use GPU + splitting
    options.RC =  RC;
    options.dBThresh = -18;
    ElementSplitting = 2;
    options.ElementSplitting = [1 ElementSplitting];
    param.bandwidth = 70;

    % GPU-distributed computation: divide into point chunks
    gpu = gpuDevice;
    disp(['GPU available ',num2str(gpu.AvailableMemory/1e6),' MB mem'])
    ReduCoef = 1.1;

    MemoryUnit = 5e6 * param.Nelements/128 * 2 * 2 * ElementSplitting;
    Nx = numel(x);  
    totalMemory = MemoryUnit * ceil(Nx/1e4);
    NW = ceil(totalMemory/(gpu.AvailableMemory/ReduCoef));
    numPointsCal = floor(gpu.AvailableMemory/ReduCoef/MemoryUnit) * 1e4;
    dim1Dist = [ones(1,NW-1)*numPointsCal Nx-(NW-1)*numPointsCal];
    x = mat2cell(x(:),dim1Dist,1);
    y = mat2cell(y(:),dim1Dist,1);
    z = mat2cell(z(:),dim1Dist,1);
    allRC = options.RC;
    RC = mat2cell(allRC(:),dim1Dist,1);

    for i = 1:NW
        disp(['--- Distributed chunk: ',num2str(i),'/',num2str(NW),' ---']);
        options.RC = RC{i};
        [~,~,RFspD{i},idx] = alg.FQpfield3(x{i},y{i},z{i},delaysTX,param,options);
    end
    RFsp = zeros(size(RFspD{1}));
    for i = 1:NW
        RFsp = RFsp + RFspD{i};
    end
    RFspectrum(idx,:) = RFsp;
end

% ------------ TIME DOMAIN RF VIA IFFT ----------------
nf = ceil(param.fs/2/param.fc*(Nf-1));
RF = ifft(conj(RFspectrum),nf,'symmetric');
RF = RF(1:floor(nf/2),:);

% Threshold very small values (denoising)
RelThresh = 1e-5; 
tmp = @(RelRF) round(0.5*(1+tanh((RelRF-RelThresh)/(RelThresh/10)))/(RelThresh/10))*(RelThresh/10);
RF = RF.*tmp(abs(RF)/(eps+max(abs(RF(:)))));

end

% ------------ FIELD NAME NORMALIZATION ----------------
function structArray = IgnoreCaseInFieldNames(structArray)
switch inputname(1)
    case 'param'
        fieldLIST = {'attenuation','baffle','bandwidth','c','fc',...
            'fnumber','focus','fs','height','kerf','movie','Nelements',...
            'passive','pitch','radius','RXangle','RXdelay',...
            'TXapodization','TXdelay','TXfreqsweep','TXnow','t0','width'};
    case 'options'
        if isstruct(structArray)
            fieldLIST = {'dBThresh','ElementSplitting',...
                'FullFrequencyDirectivity','FrequencyStep','ParPool',...
                'WaitBar'};
        else
            return
        end
end
OldFieldNames = fieldnames(structArray);
tmp = lower(OldFieldNames);
assert(length(tmp)==length(unique(tmp)),...
    ['Duplicate field names exist in ', upper(inputname(1))])
[idx,loc] = ismember(lower(fieldLIST),tmp);
idx = find(idx); loc = loc(idx);
for k = 1:length(idx)
    tmp = eval(['structArray.' OldFieldNames{loc(k)}]); %#ok
    structArray = rmfield(structArray,OldFieldNames{loc(k)});
    eval(['structArray.' fieldLIST{idx(k)} ' = tmp;']) %#ok
end
end

% ------------ EXAMPLE USAGE ----------------
function [RF,param] = RunTheExample
param.fc = 3e6;
param.bandwidth = 70;
param.width = 250e-6;
param.height = 250e-6;

pitch = 300e-6;
[xe,ye] = meshgrid(((1:32)-16.5)*pitch);
param.elements = [xe(:).'; ye(:).'];

x0 = 0; y0 = 0; z0 = 3e-2;
dels = txdelay3(x0,y0,z0,param);

N = 100;
x = 2*(rand(1,N)-0.5)*4e-3;
y = 2*(rand(1,N)-0.5)*4e-3;
z = rand(1,N)*6e-2;
RC = hypot(rand(1,N),rand(1,N));

param.fs = 10*param.fc;
[RF,param] = FQsimus3(x,y,z,RC,dels,param);

figure,scatter3(x*1e3,y*1e3,z*1e3,30,RC,'filled')
colormap(cool),hold on
scatter3(xe(:)*1e3,ye(:)*1e3,0*xe(:),3,'b','filled')
axis equal,box on,set(gca,'zdir','reverse'),zlabel('[mm]')
title([int2str(N) ' scatterers'])
end