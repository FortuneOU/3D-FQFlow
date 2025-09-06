classdef Quantification
    properties
        floderPath   % Path to result data folder
        globalParam  % Global configuration structure
    end
    methods
        % Constructor
        function obj = Quantification(globalParam)
            obj.floderPath  = fullfile('..\Result',globalParam.resultFloder);
            obj.globalParam = globalParam;
        end

        % Evaluate quantitative metrics of reconstruction
        function [] = runQuantification(obj)
           xi = obj.globalParam.ReconArea.xi;
           yi = obj.globalParam.ReconArea.yi;
           zi = obj.globalParam.ReconArea.zi;
           voxelsize = abs(xi(1) - xi(2));

           % Select phantom case and corresponding geometry transform
           switch obj.globalParam.phantomCase
                case 'Renal'
                    geometyData_path = '..\Data\25_03_14_03_10';
                    T = utils.spaceTransform([3 2 1], [11.3,6.8,9]*1e-3, [0,0,24]*1e-3, "RenalTran");
                case 'AiVsl'
                    geometyData_path = '..\Data\25_05_20_02_26';
                    T = utils.spaceTransform([1 3 2], [-1.6,8.2,34.4]*1e-3, [0,0,14]*1e-3, "AiVslTran");
                case 'Gln'
                    geometyData_path = '..\Data\25_05_20_03_02';
                    T = utils.spaceTransform([1 2 3], [0,0,0]*1e-3, [40,45,-50]*1e-3, "GlnTran");
                otherwise
                    error('Unknown case configuration');
            end

            % Load ground-truth flow geometry
            vtufilepath = fullfile(geometyData_path,'vtu.mat');
            GeometryPropertiesPath = fullfile(geometyData_path,'GeometryProperties.mat');
            load(GeometryPropertiesPath,'vtuProperties');
            [vtuStruct2, ~] = load_vessel_data(vtufilepath,vtuProperties);

            % Transform to imaging coordinates
            positions  = T.transform(vtuStruct2.points);
            velocities = T.transformVec(vtuStruct2.velocities);

            % Subsample for visualization
            dec = 1;
            sampled_positions  = positions(1:dec:end, :);
            sampled_velocities = velocities(1:dec:end, :);

            % ------ Build ground-truth voxel flow volume -----
            p_start = [xi(1),yi(1),zi(1)];
            gt_mask = false(length(xi),length(yi),length(zi));
            vel_Volume = zeros(length(xi),length(yi),length(zi));

            % Convert physical positions to voxel indices
            pts = round((sampled_positions-p_start)/voxelsize);
            Idx =  find(all(pts >= 1 & pts <= length(xi), 2));
            sampled_velocities2 = sampled_velocities(Idx,:);
            sampled_positions2  = sampled_positions(Idx,:);
            pts = pts(Idx,:);

            for i = 1:length(pts)
                vel_Volume(pts(i,1),pts(i,2),pts(i,3)) = norm(sampled_positions2(i,:));
            end
            vol1 = vel_Volume/max(vel_Volume,[],'all');

            % ------ Load reconstructed Doppler image ------
            matName = fullfile(obj.floderPath,'uPDI.mat');
            load(matName,'Img0');
            Img2 = max(Img0, -45);
            Img1 = (Img2+45)/45;
            vol2 = Img1;

            % ------ Quantitative evaluation metrics ------
            MAX_val = double(max(vol1(:)));
            MSE = mean((double(vol1(:)) - double(vol2(:))).^2);
            PSNR = 10 * log10(MAX_val^2 / MSE);
            ncc = corrcoef(vol1(:), vol2(:));
            NCC_value = ncc(1,2);
            SSIM = ssim(vol1,vol2);

            % Report results
            disp(['MSE   = ',num2str(MSE)]);
            disp(['PSNR  = ',num2str(PSNR)]);
            disp(['SSIM  = ',num2str(SSIM)]);
            disp(['NCC   = ',num2str(NCC_value)]);
        end
    end
end