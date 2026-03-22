classdef Scene
    properties
        MotionMode       % Type of motion model selected
        %ThreeDMode       % True if 3D reconstruction is used
        probeCase        % Probe type
        phantomCase      % Phantom type
        velocityFieldMode% Type of background motion field
        tissueScatters   % Tissue scatterers outside vessels
        T                % Spatial transform (geometry alignment)
    end
    methods
        function obj = Scene(globalParam)
            % Define scatterer distribution depending on probe type
            switch globalParam.probeCase
                case 'L11-5v'
                    param = getparam('L11-5v');
                    param.c = 1540;
                    [phantomPoints(:,1),phantomPoints(:,2),phantomPoints(:,3)] = ...
                        genscat([4 6 1.5]*1e-2,param); % scatterer density
                case 'Vermon'
                    param.c = 1540; param.fc = 7.81e6; param.bandwidth = 70;
                    [phantomPoints(:,1),phantomPoints(:,2),phantomPoints(:,3)] = ...
                        genscat([2.1 2.6 2.1]*1e-2,param);
                case 'Vermon-3'
                    param.c = 1540; param.fc = 3e6; param.bandwidth = 70;
                    [phantomPoints(:,1),phantomPoints(:,2),phantomPoints(:,3)] = ...
                        genscat([4.1 8.1 4.1]*1e-2,param);
                otherwise
                    error('Unknown probeCase');
            end

            % Select geometry transform based on phantom case
            switch globalParam.phantomCase
                case 'Renal'
                    flowFileFolder = '..\Data\25_03_14_03_10';
                    obj.T = utils.spaceTransform([3 2 1], ...
                        [11.3,6.8,9]*1e-3, ...
                        [0,0,24]*1e-3, "RenalTran");
                case 'AiVsl'
                    flowFileFolder = '..\Data\25_05_20_02_26';
                    obj.T = utils.spaceTransform([1 3 2], ...
                        [-1.6,8.2,34.4]*1e-3, ...
                        [0,0,14]*1e-3, "AiVslTran");
                case 'Gln'
                    flowFileFolder = '..\Data\25_05_20_03_02';
                    obj.T = utils.spaceTransform([1 2 3], ...
                        [0,0,0]*1e-3, ...
                        [40,45,-50]*1e-3, "GlnTran");
                case 'T1Trace'
                    flowFileFolder = '..\Data\T1Trace';
                    obj.T = utils.spaceTransform([1 2 3], ...
                        [0,0,75]*1e-3, ...
                        [0,-pi/6,0], ...
                        [0,0,25]*1e-3, "T1Trace");
                case 'T2Trace'
                    flowFileFolder = '..\Data\T2Trace';
                    obj.T = utils.spaceTransform([1 3 2], ...
                        [0,0,0]*1e-3, ...
                        [0,0,0], ...
                        [7,60,55]*1e-3, "T2Trace");
                case 'T2Trace_1'
                    flowFileFolder = '..\Data\T2Trace';
                    obj.T = utils.spaceTransform([1 3 2], ...
                        [0,0,0]*1e-3, ...
                        [0,0,0], ...
                        [7,50,42]*1e-3, "T2Trace_1");
                case 'T3Trace'
                    flowFileFolder = '..\Data\T3Trace';
                    obj.T = utils.spaceTransform([1 3 2], ...
                        [0,0,0]*1e-3, ...
                        [0,0,0], ...
                        [7,50,42]*1e-3, "T2Trace_1");
                 case 'F2Trace'
                   flowFileFolder = '..\Data\F2Trace';
                    obj.T = utils.spaceTransform([1 2 3], ...
                        [0,0,0]*1e-3, ...
                        [0,0,0], ...
                        [0,25,20]*1e-3, "F2Trace");
                 case 'F20Trace'
                   flowFileFolder = '..\Data\F20Trace';
                    obj.T = utils.spaceTransform([1 2 3], ...
                        [0,0,0]*1e-3, ...
                        [0,0,0], ...
                        [0,25,20]*1e-3, 0.1,"F20Trace");
                otherwise
                    error('Unknown phantomCase');
            end

            % Load geometry information and mask out vessels
            vtufilepath = fullfile(flowFileFolder,'vtu.mat');
            GeometryPropertiesPath = fullfile(flowFileFolder,'GeometryProperties.mat');
            load(GeometryPropertiesPath,'vtuProperties');
            [~, Grid] = alg.load_vessel_data_FQ(vtufilepath,vtuProperties);

            % Transform scatterer positions to vessel grid coordinates
            phantomPointsInGrid = obj.T.invTransform(phantomPoints);
            switch globalParam.phantomCase
                case  'T1Trace'

                    [~,outsideGrid] = alg.get_vtu_indices4vtu(phantomPointsInGrid,Grid);
                    obj.tissueScatters = phantomPoints(outsideGrid,:);
                case  'T2Trace'
                    [~,outsideGrid] = alg.get_vtu_indices4vtu(phantomPointsInGrid,Grid);
                    obj.tissueScatters = phantomPoints(outsideGrid,:);

                otherwise
                    inVesselFlag = alg.get_vtu_indices(phantomPointsInGrid,Grid);
                    obj.tissueScatters = phantomPoints((inVesselFlag == 0),:);
            end
            % Keep only tissue scatterers outside vessel lumen


            % Store setup params
            obj.MotionMode       = globalParam.MotionMode;
            obj.velocityFieldMode= globalParam.velocityFieldMode;
            obj.probeCase        = globalParam.probeCase;
            obj.phantomCase      = globalParam.phantomCase;
        end
    end
end
