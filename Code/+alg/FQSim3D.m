function [] = FQSim3D(scatterPos,scatterRC,outputName,param)
% Use MUST toolbox to simulate ultrasound data (3D case)

% Define point scatterers
xs = scatterPos(:,1); % x-coordinates [m]
ys = scatterPos(:,2); % y-coordinates [m]
zs = scatterPos(:,3); % z-coordinates [m]
RC = scatterRC;       % Reflection coefficients

% Load transmit delay laws
txdel = param.txdel;
Na = param.Na; % Number of transmits

% Generate RF signals
RF = cell(Na,1);       % RF signals for each transmit event
param.fs = 4*param.fc; % Sampling frequency (Hz)
option.WaitBar = true; % Show waitbar
param.c = 1540;
for k = 1:Na
    disp(['Simulating TX event: ',num2str(k)]);
    tic;
    % 3D simulation (using modified GPU-accelerated version of simus3)
    RF{k} = alg.FQsimus3(xs,ys,zs,RC,txdel{k},param,option);
    elapsedTime = toc;
    disp(['This TX required ', num2str(elapsedTime/60, '%.2f'), ' minutes']);
end

% Save results
folderPath = fileparts(outputName);
if ~isfolder(folderPath)
    mkdir(folderPath);
end
save(outputName, 'RF','-v7.3'); % Use v7.3 format for large datasets
end