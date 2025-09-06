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
                        [0,0,24]*1e-3, "RenalTran");
                case 'AiVsl'
                    obj.flowFileFolder = '..\Data\25_05_20_02_26';
                    obj.T = utils.spaceTransform([1 3 2], ...
                        [-1.6,8.2,34.4]*1e-3, ...
                        [0,0,14]*1e-3, "AiVslTran");
                case 'Gln'
                    obj.flowFileFolder = '..\Data\25_05_20_03_02';
                    obj.T = utils.spaceTransform([1 2 3], ...
                        [0,0,0]*1e-3, ...
                        [40,45,-50]*1e-3, "GlnTran");
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
                tempPos = Frame.Pulse1.Points(1:100:end,:);  % scatterer positions
                tempRc  = Frame.Pulse1.Radius(1:100:end); % scatterer radii
                obj.flowScatters{k}.pos = obj.T.transform(tempPos);                
                % Reflection coefficient scaling
                obj.flowScatters{k}.rc  = (tempRc./8e-7).^6 * obj.RCcoef;
            end


        end
    end
end