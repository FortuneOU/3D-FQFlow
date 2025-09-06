function [] = FQBmFr2D(inputName,outputName,param,Area)
load(inputName);
param.passive = true;
Na = param.Na;
txdel = param.txdel2d; % Transmit delay laws for 2D configuration % Demodulation of RF signals IQ = cell(Na,1); % Container for baseband I/Q signals
for k = 1:Na
    IQ{k} = rf2iq(RF{k},param.fs,param.fc);
end % Delay-and-Sum Beamforming
param.fnumber = 1.7;
param.c = 1540;
param.t0 = 0;
%xi = squeeze(Area.xi(:,1,:)); zi = squeeze(Area.zi(:,1,:));
[xi,zi] = meshgrid(linspace(-2e-2,2e-2,200),linspace(0.5e-2,5e-2,225));

bIQ = zeros(length(zi),width(zi),Na); % Stores I/Q images for each transmit event
h = waitbar(0,'');
for k = 1:Na
    waitbar(k/Na,h,['Beamforming progress: I/Q dataset #' int2str(k) ' of ',num2str(Na)])
    bIQ(:,:,k) = das(IQ{k},xi,zi,txdel{k},param);
end
close(h)
ReconIQ = mean(bIQ,3); % Final compound image by averaging all transmissions % Save reconstruction results
folderPath = fileparts(outputName);
if ~isfolder(folderPath)
    mkdir(folderPath);
end
save(outputName, 'ReconIQ','bIQ');
