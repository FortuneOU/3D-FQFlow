classdef ImageProcessor
    properties
        floderPath   % Path to store IQ data files
        ReconArea    % Reconstruction area definition
        globalParam  % Global parameter configuration
    end
    methodsclassdef ImageProcessor
    properties
        floderPath   % Path to store IQ data files
        ReconArea    % Reconstruction area definition
        globalParam  % Global parameter configuration
        param        % us parameeter

    end
    methods
        % Constructor
        function obj = ImageProcessor(globalParam,usParam)
            obj.floderPath  = fullfile('..\Result',globalParam.resultFloder,'IQ');
            obj.ReconArea   = globalParam.ReconArea;
            obj.globalParam = globalParam;
            obj.param.fs = usParam.fs;
            obj.param.fc = usParam.fc;
            obj.param.c = 1540;
            obj.param.bandwidth = usParam.bandwidth;
            obj.param.width = usParam.width;
            obj.param.height = usParam.height;
            obj.param.radius = usParam.radius;
            obj.param.Nelements = usParam.Nelements;
            obj.param.pitch = usParam.pitch;
            obj.param.focus = usParam.focus;
            obj.param.PRF = usParam.PRF;
            obj.param.elements = usParam.elements;
            obj.param.fnumber = usParam.fnumber;
            obj.param.t0 = usParam.t0;
            obj.param.passive = usParam.passive;
            obj.param.Na = usParam.Na;
            obj.param.txdel2d = usParam.txdel2d;


        end

        % 2D image post-processing with SVD clutter filter
        function [] = process2DImages4BanduPDI(obj)
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
            for i = 1:200
                filePath = fullfile(obj.floderPath, sortedFiles(i).name);
                fileData = load(filePath);
                IQData(:,:,i) = fileData.ReconIQ;
            end

            % Perform SVD clutter filtering
            packetSize = numFiles;
            SVD_input_2d0 = reshape(IQData,[],packetSize);
            [U, S, V] = svd(SVD_input_2d0,'econ');

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
            x_image = squeeze(obj.ReconArea.xi(1,:,1));

            y_image = squeeze(obj.ReconArea.zi(1,1,:));
            Img0 = db(I./max(I,[],"all"));

            I_uPDI_80db = max(Img0,-80);

            % uPDI image visualization
            figure;
            imagesc(x_image,y_image,Img0);
            caxis([-60 0]); axis image; colormap hot;
            title('Ultrasound Power Doppler (uPDI) Image');

            % B-mode reference image (from first frame)
            B_mode_1 = abs(IQData(:,:,1));
            B_mode   = db(B_mode_1./max(B_mode_1,[],"all"));

            I_B_120db = max(B_mode,-120);

            figure;
            imagesc(x_image,y_image,B_mode);
            caxis([-45 0]); axis image; colormap gray;
            title('B-mode Image');

            % Save processed Doppler image
            % matName = fullfile('..\Result',obj.globalParam.resultFloder,'processedImage.mat');
            % save(matName,'I_uPDI_80db','I_B_120db');
        end

        % 3D image post-processing with SVD clutter filter
        function [] = process3DImages(obj)
            HPF = 2; % Threshold rank index for SVD
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
            %
            % SVD filtering
            packetSize = numFiles;
            IQData2 = permute(IQData,[2,1,3,4]);
            SVD_input_2d0 = reshape(IQData2,[],packetSize);
            SVD_input_2d0 = single(SVD_input_2d0);
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

            B_mode_1 = abs(IQData(:,:,:,1));
            B_mode   = db(B_mode_1./max(B_mode_1,[],"all"));

            I_B_120db = max(B_mode,-120);
            I_uPDI_80db = max(Img0,-80);

            % Save volumetric Doppler image
            matName = fullfile('..\Result',obj.globalParam.resultFloder,'processedImage.mat');
            save(matName,'I_uPDI_80db','I_B_120db');
        end

        % 2D image post-processing FOR COLOR DOPPLER
        function [] = process2DImages4C(obj)

            saveSeqFlag = true;

            fileList = dir(fullfile(obj.floderPath, 'IQ_*.mat'));
            packetsize_C = 40*2;
            obj.param.PRF = obj.param.PRF/2;
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
            IQData = zeros(m, n, packetsize_C);
            start_num = 0;
            if saveSeqFlag == true
                loopNum = floor((numFiles - start_num + 1)/packetsize_C);

                for k = 1:loopNum
                    disp(['Calc ',num2str(k),'th sets'])
                    % Sequential load of IQ frames
                    
                    for i = start_num+1:2:start_num+packetsize_C
                        filePath = fullfile(obj.floderPath, sortedFiles(i).name);
                        fileData = load(filePath);
                        IQData(:,:,ceil((i-start_num)/2)) = fileData.ReconIQ;
                    end
                    start_num = start_num+packetsize_C;
                    Vd(:,:,k) = iq2doppler(IQData,obj.param,[7 7]);
                end
                % Save processed Doppler image
                matName = fullfile('..\Result',obj.globalParam.resultFloder,'CmodeImages.mat');
                save(matName,'Vd');

            else

                % Sequential load of IQ frames
                for i = start_num+1:start_num+packetsize_C
                    filePath = fullfile(obj.floderPath, sortedFiles(i).name);
                    fileData = load(filePath);
                    IQData(:,:,i-start_num) = fileData.ReconIQ;
                end

                Vd = iq2doppler(IQData,obj.param,[7 7]);

                x_image = squeeze(obj.ReconArea.xi(1,:,1));
                y_image = squeeze(obj.ReconArea.zi(1,1,:));
                figure;
                imagesc(x_image*100,y_image*100,Vd);

                colormap dopplermap
                colorbar
                caxis([-1 1]*max(abs(Vd(:))));

                title('Color Doppler map ')
                ylabel('[cm]')
                axis equal ij tight
                set(gca,'XColor','none','box','off')
            end
            % Save processed Doppler image
            % matName = fullfile('..\Result',obj.globalParam.resultFloder,'processedImage.mat');
            % save(matName,'I_uPDI_80db','I_B_120db');
        end

        % 2D image post-processing FOR COLOR DOPPLER
        function [] = process2DImages4PW(obj)

            fileList = dir(fullfile(obj.floderPath, 'IQ_*.mat'));
            packetsize_PW = 4000;
            saveSeqFlag = true;
            obj.param.PRF = obj.param.PRF/2;

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
            IQData = zeros(m, n, packetsize_PW/2);

            start_num = 0;

            % Sequential load of IQ frames
            for i = start_num+1:2:2:start_num+packetsize_PW
                filePath = fullfile(obj.floderPath, sortedFiles(i).name);
                fileData = load(filePath);
                IQData(:,:,ceil((i-start_num)/2)) = fileData.ReconIQ;
            end


            %IQ_gate = squeeze(mean(mean(IQData(209:211,106:108,:),1),2));
            %IQ_gate = squeeze(mean(mean(IQData(125:127,44:46,:),1),2));
            %IQ_gate = squeeze(mean(mean(IQData(80:82,73:75,:),1),2));
            IQ_gate = squeeze(mean(mean(IQData(130:132,138:140,:),1),2));

            %% ====== 多普勒频谱计算 ======
            N = length(IQ_gate);             % 脉冲数量
            win_len = 128;                    % STFT窗长度(多普勒窗，例如64脉冲)
            overlap = 100;                    % 窗重叠
            window = hanning(win_len);
            theta = 0;


            % 参数
            % fs = obj.param.fs;
            fcn = 0.1;

            % 设计高通FIR滤波器
            N = 32;
            b = fir1(N, fcn, 'high', hamming(N+1));
            IQ_fir = filtfilt(b, 1, IQ_gate);

            % 使用短时傅里叶变换计算时间-频率谱
            [S, f_d, t_s] = spectrogram(IQ_fir, window, overlap, win_len, obj.param.PRF, 'centered');

            % 计算功率谱
            S_power = abs(S).^2;
            % S_power(floor(win_len/2)-floor(win_len/20):floor(win_len/2)+floor(win_len/20),:) = 0;

            % 将频率轴转换为速度
            velocity_axis = (obj.param.c * f_d) / (4*obj.param.fc*cos(theta));   % m/s

            %% ====== 绘图 ======
            figure;
            imagesc(t_s, velocity_axis, 10*log10(S_power./max(S_power(:))));
            axis xy;
            colormap gray;
            clim([-15 0])
            xlabel('Time (s)');
            ylabel('Velocity (m/s)');
            title('PW Doppler Spectrum');
            colorbar;
        end


    end
end
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
