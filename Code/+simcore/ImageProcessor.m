classdef ImageProcessor
    properties
        floderPath   % Path to store IQ data files
        ReconArea    % Reconstruction area definition
        globalParam  % Global parameter configuration
    end
    methods
        % Constructor
        function obj = ImageProcessor(globalParam)
            obj.floderPath  = fullfile('..\Result',globalParam.resultFloder,'IQ');
            obj.ReconArea   = globalParam.ReconArea;
            obj.globalParam = globalParam;
        end

        % 2D image post-processing with SVD clutter filter
        function [] = process2DImages(obj)
            HPF = 20; % Threshold for SVD-based clutter rejection
            fileList = dir(fullfile(obj.floderPath, 'IQ_*.mat'));

            % Extract numeric indices from filenames
            fileNames = {fileList.name};
            numFiles  = numel(fileNames);
            fileNumbers = zeros(1, numFiles);
            for k = 1:numFiles
                match = regexp(fileNames{k}, 'IQ_(\d+)\.mat', 'tokens');
                fileNumbers(k) = str2double(match{1}{1});
            end
            [~, sortedIdx] = sort(fileNumbers);
            sortedFiles = fileList(sortedIdx);

            % Load first file to preallocate
            sampleData = load(fullfile(obj.floderPath, sortedFiles(1).name));
            [m, n] = size(sampleData.ReconIQ);
            IQData = zeros(m, n, numFiles);

            % Sequential load of IQ frames
            for i = 1:numFiles
                filePath = fullfile(obj.floderPath, sortedFiles(i).name);
                fileData = load(filePath);
                IQData(:,:,i) = fileData.ReconIQ;
            end

            % Perform SVD clutter filtering
            packetSize = numFiles;
            SVD_input_2d0 = reshape(IQData,[],packetSize);
            [U, S, V] = svd(SVD_input_2d0,'econ');
            diag(S)
            % Display correlation matrix for tissue suppression visualization
            figure;
            TSMatrix = corrcoef(abs(U));
            imagesc(TSMatrix),colormap(jet),colorbar;
            caxis([0 1]);

            % Construct filter: keep high-rank singular values (blood)
            v = zeros(packetSize,1);
            v(HPF:packetSize) = 1;
            vv = diag(v);
            S = S.*vv;

            % Apply filter and reshape
            SVD_output_2d0 = U*S*V';
            Result0 = reshape(SVD_output_2d0,size(IQData));
            I = mean(abs(Result0).^2,3);

            % Normalize image in dB
            x_image = squeeze(obj.ReconArea.xi(:,1,:));
            y_image = squeeze(obj.ReconArea.zi(:,1,:));
            Img0 = db(I./max(I,[],"all"));

            % uPDI image visualization
            figure;
            imagesc(x_image,y_image,Img0);
            caxis([-60 0]); axis image; colormap hot;
            title('Ultrasound Power Doppler (uPDI) Image');

            % B-mode reference image (from first frame)
            B_mode_1 = abs(IQData(:,:,1));
            B_mode   = db(B_mode_1./max(B_mode_1,[],"all"));
            figure;
            imagesc(x_image,y_image,B_mode);
            caxis([-45 0]); axis image; colormap gray;
            title('B-mode Image');

            % Save processed Doppler image
            matName = fullfile('..\Result',obj.globalParam.resultFloder,'uPDI.mat');
            save(matName,Img0);
        end

        % 3D image post-processing with SVD clutter filter
        function [] = process3DImages(obj)
            HPF = 20; % Threshold rank index for SVD
            fileList = dir(fullfile(obj.floderPath, 'IQ_*.mat'));

            % Extract frame indices
            fileNames = {fileList.name};
            numFiles = numel(fileNames);
            fileNumbers = zeros(1, numFiles);
            for k = 1:numFiles
                match = regexp(fileNames{k}, 'IQ_(\d+)\.mat', 'tokens');
                fileNumbers(k) = str2double(match{1}{1});
            end
            [~, sortedIdx] = sort(fileNumbers);
            sortedFiles = fileList(sortedIdx);

            % Load one volume to preallocate
            sampleData = load(fullfile(obj.floderPath, sortedFiles(1).name));
            [m, n, q] = size(sampleData.ReconIQ);
            IQData = zeros(m, n, q,numFiles);

            % Load all IQ volumes
            for i = 1:numFiles
                filePath = fullfile(obj.floderPath, sortedFiles(i).name);
                fileData = load(filePath);
                IQData(:,:,:,i) = fileData.ReconIQ;
            end

            % SVD filtering
            packetSize = numFiles;
            IQData2 = permute(IQData,[2,1,3,4]);
            SVD_input_2d0 = reshape(IQData2,[],packetSize);
            [U, S, V] = svd(SVD_input_2d0,'econ');

            % Correlation display
            figure;
            TSMatrix = corrcoef(abs(U));
            imagesc(TSMatrix),colormap(jet),colorbar;

            % Reject low-rank components (clutter)
            v = zeros(packetSize,1);
            v(HPF:packetSize) = 1;
            vv = diag(v);
            S = S.*vv;

            % Apply filter and reshape result
            SVD_output_2d0=U*S*V';
            Result0=reshape(SVD_output_2d0,size(IQData));
            I = mean(abs(Result0).^2,4);

            % Normalize in dB
            x_image = obj.ReconArea.xi;
            y_image = obj.ReconArea.yi;
            z_image = obj.ReconArea.zi;
            Img0 = db(I./max(I,[],"all"));

            % Save volumetric Doppler image
            matName = fullfile('..\Result',obj.globalParam.resultFloder,'uPDI.mat');
            save(matName,'Img0');
        end
    end
end
