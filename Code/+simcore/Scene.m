classdef Scene
    properties
        MotionMode       % Type of motion model selected
        ThreeDMode       % True if 3D reconstruction is used
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
                        genscat([4 5 0.2]*1e-2,param); % scatterer density
                case 'Vermon'
                    param.c = 1540; param.fc = 7.81e6; param.bandwidth = 70;
                    [phantomPoints(:,1),phantomPoints(:,2),phantomPoints(:,3)] = ...
                        genscat([2.1 2.6 2.1]*1e-2,param);
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
                otherwise
                    error('Unknown phantomCase');
            end

            % Load geometry information and mask out vessels
            vtufilepath = fullfile(flowFileFolder,'vtu.mat');
            GeometryPropertiesPath = fullfile(flowFileFolder,'GeometryProperties.mat');
            load(GeometryPropertiesPath,'vtuProperties');
            [~, Grid] = alg.load_vessel_data(vtufilepath,vtuProperties);

            % Transform scatterer positions to vessel grid coordinates
            phantomPointsInGrid = obj.T.invTransform(phantomPoints);
            inVesselFlag = alg.get_vtu_indices(phantomPointsInGrid,Grid);

            % Keep only tissue scatterers outside vessel lumen
            obj.tissueScatters = phantomPoints((inVesselFlag == 0),:);

            % Store setup params
            obj.MotionMode       = globalParam.MotionMode;
            obj.velocityFieldMode= globalParam.velocityFieldMode;
            obj.probeCase        = globalParam.probeCase;
            obj.phantomCase      = globalParam.phantomCase;
        end
    end
end