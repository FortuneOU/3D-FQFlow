classdef Scatterers
    properties
        probeCase       % Probe type (e.g., L11-5v, Vermon)
        phantomCase     % Phantom type (e.g., Renal, AiVsl, Gln)
        flowScatters    % Flow-related scatterers (cell per frame)
        tissueScatters  % Tissue scatterers
        T               % Space transform object
        flowFileFolder  % Folder containing CFD-derived flow scatterers
        RCcoef          % Reflection coefficient scaling factor
    end
    methods
        function obj = Scatterers(globalParam)
            % Select geometry folder and coordinate transform based on phantomCase
            switch globalParam.phantomCase
                case 'Renal'
                    obj.flowFileFolder = '..\Data\25_03_14_03_10';
                    obj.T = utils.spaceTransform([3 2 1], ...
                        [11.3,6.8,9]*1e-3, ...
                        [0,0,0], ...
                        [0,0,24]*1e-3, "RenalTran");
                case 'AiVsl'
                    obj.flowFileFolder = '..\Data\25_05_20_02_26';
                    obj.T = utils.spaceTransform([1 3 2], ...
                        [-1.6,8.2,34.4]*1e-3, ...
                        [0,0,0], ...
                        [0,0,14]*1e-3, "AiVslTran");
                case 'Gln'
                    obj.flowFileFolder = '..\Data\25_05_20_03_02';
                    obj.T = utils.spaceTransform([1 2 3], ...
                        [0,0,0]*1e-3, ...
                        [0,0,0], ...
                        [40,45,-50]*1e-3, "GlnTran");

                case 'T1Trace'
                    obj.flowFileFolder = '..\Data\T1Trace';
                    obj.T = utils.spaceTransform([1 2 3], ...
                        [0,0,75]*1e-3, ...
                        [0,-pi/6,0], ...
                        [0,0,25]*1e-3, "T1Trace");
                case 'T2Trace'
                    obj.flowFileFolder = '..\Data\T2Trace';
                    obj.T = utils.spaceTransform([1 3 2], ...
                        [0,0,0]*1e-3, ...
                        [0,0,0], ...
                        [7,60,55]*1e-3, "T2Trace");
                case 'T2Trace_1'
                    obj.flowFileFolder = '..\Data\T2Trace';
                    obj.T = utils.spaceTransform([1 3 2], ...
                        [0,0,0]*1e-3, ...
                        [0,0,0], ...
                        [7,50,42]*1e-3, "T2Trace_1");
                case 'T3Trace'
                    obj.flowFileFolder = '..\Data\T3Trace';
                    obj.T = utils.spaceTransform([1 3 2], ...
                        [0,0,0]*1e-3, ...
                        [0,0,0], ...
                        [7,50,42]*1e-3, "T2Trace_1");
                 case 'F2Trace'
                    obj.flowFileFolder = '..\Data\F2Trace';
                    obj.T = utils.spaceTransform([1 2 3], ...
                        [0,0,0]*1e-3, ...
                        [0,0,0], ...
                        [0,25,20]*1e-3, "F2Trace");
                 case 'F20Trace'
                    obj.flowFileFolder = '..\Data\F20Trace';
                    obj.T = utils.spaceTransform([1 2 3], ...
                        [0,0,0]*1e-3, ...
                        [0,0,0], ...
                        [0,0,10]*1e-3, 0.1,"F20Trace");
                otherwise
                    error('Unknown phantomCase');
            end
            obj.probeCase = globalParam.probeCase;
            obj.RCcoef    = globalParam.flowRCcoef;
        end

        % Load flow scatterer positions and convert into imaging coordinates
        function obj = loadFlowScatters(obj)
            fileList = dir(fullfile(obj.flowFileFolder, 'Frame_*.mat'));
            for k = 1:length(fileList)
                disp(['Loading frame #',num2str(k),' flow scatterers...']);
                io = utils.IO(fileList(k).name, obj.flowFileFolder);
                io = io.readData('Frame');
                Frame = io.data;
                tempPos = Frame.Pulse1.Points(1:10000:end,:);  % scatterer positions
                tempRc  = Frame.Pulse1.Radius(1:10000:end); % scatterer radii
                obj.flowScatters{k}.pos = obj.T.transform(tempPos);
                % Reflection coefficient scaling
                obj.flowScatters{k}.rc  = (tempRc./8e-7).^6 * obj.RCcoef;
            end
        end

        % function obj = loadFlowScattersFromMat(obj)
        %     fileName = dir(fullfile(obj.flowFileFolder, 'particle_trajectories.mat'));
        %
        %     disp('Loading flow scatterers... from Mat');
        %     io = utils.IO('particle_trajectories_10.mat', obj.flowFileFolder);
        %     io = io.readData('scatters');
        %     scatters = io.data;
        %     for i = 1:200
        %         tempPos = scatters{i+800}.pos(1:1:end,:)/100;  % scatterer positions
        %         obj.flowScatters{i}.pos = obj.T.transform(tempPos);
        %         % Reflection coefficient scaling
        %         obj.flowScatters{i}.rc  = ones(length(tempPos),1) * obj.RCcoef;
        %     end
        % end

        function obj = loadFlowScattersFromMat4Seq(obj)
            disp('Loading flow scatterers... from multiple Mat files');

            scatterCount = 0; % 计数器用于存放到 obj.flowScatters


            for fileIdx = 1:20

                matFileName = sprintf('particle_trajectories_cycle03_part%02d.mat', fileIdx);
                fprintf('Reading %s\n', matFileName);


                io = utils.IO(matFileName, obj.flowFileFolder);
                io = io.readData('scatters');
                scatters = io.data;


                %scatterRanges = [1:200, 501:700];
                scatterRanges = [1:25];
                
                for sIdx = scatterRanges
                    scatterCount = scatterCount + 1;

                    tempPos = scatters{sIdx}.pos(1:end, :) / 100;   % scatterer positions
                    obj.flowScatters{scatterCount}.pos = obj.T.transform(tempPos);

                    % Reflection coefficient scaling
                    obj.flowScatters{scatterCount}.rc  = ones(length(tempPos), 1) * obj.RCcoef;
                end
            end

        end

    end
end
