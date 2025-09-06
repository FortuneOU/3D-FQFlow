function [] = FQBmFr3D(inputFloder,outputFloder,param,Area)
% Fast 3D beamforming with multiple frames
% Compute transmit delay matrices only once, reuse for all frames

fileList = dir(fullfile(inputFloder, 'RF_*.mat'));
load(fullfile(fileList(1).folder, fileList(1).name));
isIQ = ~isreal(RF{1}); % Check if raw input already is I/Q signals
if ~isIQ
    isIQ = 1;
end
if ~isfield(param,'radius')
    param.radius = Inf; % Default: linear matrix array
else
    assert(isinf(param.radius),...
        'DASMTX3 does not handle convex arrays.')
end

%% Define interpolation scheme
method = 'linear';
tmp = strcmpi(method,...
    {'nearest','linear','quadratic','lanczos3','5points','lanczos5'});
if ~any(tmp)
    error(['METHOD must be ''nearest'', ''linear'', ''quadratic'',',...
        ' ''Lanczos3'', ''5points'' or ''Lanczos5''.'])
end
Npoints = find(tmp);

xi = Area.xi;
zi = Area.zi;
yi = Area.yi;

% ---- Memory-aware chunking ----
% The DAS matrix can become too large for memory. Split into chunks.
if ispc
    mem = memory;
    MPAB = mem.MaxPossibleArrayBytes;
end
NoE = size(RF{1},2); % Number of transducer elements
bytes = 16*NoE*Npoints*numel(xi);
factor = 5; % Adjustment factor to balance loops vs. RAM
if ispc
    Nchunks = ceil(factor*bytes*param.Na/MPAB);
else
    Nchunks = 10; % Static value for Linux/Mac
end
% -------------------------------

idx = round(linspace(1,numel(xi)+1,Nchunks+1));
delaysTX = param.txdel{1};
[nl2,nc,~] = size(RF{1});
nl = floor(nl2*1.5);
Na = param.Na;
clear RF;

for chunk_i = 1:Nchunks
    xi2 = xi(idx(chunk_i):idx(chunk_i+1)-1);
    yi2 = yi(idx(chunk_i):idx(chunk_i+1)-1);
    zi2 = zi(idx(chunk_i):idx(chunk_i+1)-1);
    for tx_k = 1:Na
        disp(['Processing chunk:',num2str(chunk_i),'/',num2str(Nchunks),...
            '    TX:',num2str(tx_k),'/',num2str(Na)])
        var1{1} = param.txdel{tx_k};
        var1{2} = param;
        var1{3} = method;
        [M{tx_k},~] = FQdasmtx3((~isIQ*1+isIQ*1i)*[nl nc],...
            xi2,yi2,zi2,var1{:});
    end
    for frame_j = 1:length(fileList)
        load(fullfile(fileList(frame_j).folder, fileList(frame_j).name)) ;
        for tx_j = 1:Na
            b = rf2iq(RF{tx_j},param.fs,param.fc);
            padRows = nl - size(b,1);
            SIG = [b; zeros(padRows, size(b, 2), 'like', b)];
            SIG2 = SIG(:);
            IQchunk(:,frame_j,tx_j) = M{tx_j}*SIG2;
        end
    end
    IQchunk_sumTX = mean(IQchunk,3);
    outputName = fullfile(outputFloder,['Chunk_',num2str(chunk_i,'%03d'),'.mat']);
    if ~exist(outputFloder, 'dir')
        mkdir(outputFloder);
    end
    save(outputName,'IQchunk_sumTX')
    clear IQchunk;
end

% ---- Final reconstruction by merging chunks ----
prefix = 'Chunk_';
columnCount = length(fileList);
filePattern = fullfile(outputFloder, [prefix, '*.mat']);
matFiles = dir(filePattern);
if isempty(matFiles)
    error('No chunk files were generated with prefix %s', prefix);
end

for colIdx = 1:columnCount
    dataList = [];
    for k = 1:length(matFiles)
        baseFileName = matFiles(k).name;
        fullFileName = fullfile(outputFloder, baseFileName);
        fprintf('Reading %s\n', fullFileName);
        matData = load(fullFileName);
        if isfield(matData, 'IQchunk_sumTX')
            data = matData.IQchunk_sumTX(:, colIdx);
            dataList = [dataList; data];  
        else
            warning('Variable IQchunk_sumTX not found in file: %s', fullFileName);
        end
    end
    if length(dataList) ~= length(xi)*length(yi)*length(zi)
        error('Mismatch between reconstructed vector length and image grid dimensions.')
    end
    ReconIQ = reshape(dataList, [length(zi), length(xi), length(yi)]);
    outputFileName = sprintf('IQ_%03d.mat', colIdx);
    outputFullFileName = fullfile(outputFloder, outputFileName);
    save(outputFullFileName, 'ReconIQ');
    fprintf('Saved final image: %s\n', outputFullFileName);
end
end