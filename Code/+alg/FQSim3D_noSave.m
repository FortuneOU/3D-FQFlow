function [RFout] = FQSim3D_noSave(scatterPos,scatterRC,outputName,param)
% Variant of FQSim3D without saving results directly to disk
% Returns RF signals in memory for further processing

% Define scatterers
xs = scatterPos(:,1); % x-coordinates [m]
ys = scatterPos(:,2); % y-coordinates [m]
zs = scatterPos(:,3); % z-coordinates [m]
RC = scatterRC;       % Reflection coefficients

% Load transmit delay laws
txdel = param.txdel;
Na = param.Na;

% Compute RF signals
RF = cell(Na,1);       % Container for RF waveforms
param.fs = 4*param.fc; % Sampling frequency [Hz]
option.WaitBar = true; % With waitbar
param.c = 1540;
for k = 1:Na
    disp(['Simulating TX event: ',num2str(k)]);
    tic;
    % Run 3D simulation
    RF{k} = alg.FQsimus3(xs,ys,zs,RC,txdel{k},param,option);
    elapsedTime = toc;
    disp(['This TX required ', num2str(elapsedTime/60, '%.2f'), ' minutes']);
end

% Do not save -> instead return results directly
RFout = RF;
end