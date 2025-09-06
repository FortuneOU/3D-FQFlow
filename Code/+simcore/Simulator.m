classdef Simulator
    properties
        SimuTissueMotion % Flag if tissue motion is considered
        vesselMotion     % If vessels move together with tissue
        MotionMode       % Motion model selection
        SimThreeDMode    % Simulation mode (2D/3D)
        ReconThreeDMode  % Reconstuction Mode (2D/3D)
        velocityFieldMode% Background velocity field type
        totalScatterers  % Combined scatterers
        usParam          % Ultrasound transducer parameters
        globalParam      % Global experiment parameters
        flowScatters     % Flow scatterer positions/rc per frame
        scene            % Scene object (tissue background definition)
    end
    methods
        % Constructor
        function obj = Simulator(flowScatters,scene,globalParam,usParam)
            obj.MotionMode       = globalParam.MotionMode;
            obj.SimThreeDMode    = globalParam.SimThreeDMode;
            obj.ReconThreeDMode  = globalParam.ReconThreeDMode; 
            obj.velocityFieldMode= globalParam.velocityFieldMode;
            obj.globalParam      = globalParam;
            obj.usParam          = usParam;
            obj.flowScatters     = flowScatters;
            obj.scene            = scene;

            % Simulate tissue motion trajectory and save to disk
            outputFolder = '..\Motion\MovedPhantoms';
            numPoints    = length(scene.tissueScatters);
            dt           = 1/usParam.PRF; % time step
            numFrames    = length(flowScatters);
            points       = scene.tissueScatters;

            if isempty(dir(outputFolder))
                mkdir(outputFolder);
                switch globalParam.velocityFieldMode
                    case 1 % Uniform field
                        velocity_field = @(pts) [2e-3*ones(numPoints,1), zeros(numPoints,1), 4e-3*ones(numPoints,1)];
                        for step = 1:numFrames
                            velocities = velocity_field(points);
                            points = points + velocities * dt;
                            TissuePoints = points;
                            save(fullfile(outputFolder,['MovedPhy_',num2str(step),'.mat']), 'TissuePoints');
                        end
                    case 2 % Rotational field
                        velocity_field = @(pts) [pts(:,3), zeros(numPoints,1), -pts(:,1)];
                        for step = 1:numFrames
                            velocities = velocity_field(points);
                            points = points + velocities * dt;
                            TissuePoints = points;
                            save(fullfile(outputFolder,['MovedPhy_',num2str(step),'.mat']), 'TissuePoints');
                        end
                    case 3 % Realistic field from video motion estimation
                        % Use optical flow between frames of reference movie
                        v = VideoReader('..\Motion\Kidney\202503152031937.avi');
                        try
                            frames = read(v);
                        catch
                            frames = [];
                            while hasFrame(v)
                                frames = cat(4, frames, readFrame(v));
                            end
                        end
                        frameA = rgb2gray(squeeze(frames(62:469,351:611,:,301)));
                        frameB = rgb2gray(squeeze(frames(62:469,351:611,:,325)));
                        imgA = imresize(frameA, [600, 400]);
                        imgB = imresize(frameB, [600, 400]);

                        opticFlow = opticalFlowHS('Smoothness', 0.01, 'MaxIteration', 500);
                        flow = estimateFlow(opticFlow, imgA);
                        flow = estimateFlow(opticFlow, imgB);
                        Vx = medfilt2(flow.Vx, [5 5]);
                        Vz = medfilt2(flow.Vy, [5 5]);

                        velocityField(:,:,1) = Vx*1e-3;
                        velocityField(:,:,2) = Vz*1e-3;
                        save ..\Motion\VelocityField\RenalMotion.mat velocityField

                        for step = 1:numFrames
                            velocities = velocity_field2(velocityField,points);
                            points = points + velocities * dt;
                            TissuePoints = points;
                            save(fullfile(outputFolder,['MovedPhy_',num2str(step),'.mat']), 'TissuePoints');
                        end
                    otherwise
                        obj.MotionMode = 0;  % Default: no motion
                end
            end
        end

        % Main simulation
        function [] = runSimulation(obj)
            numFrames = length(obj.flowScatters);
            outputFolder = '..\Motion\MovedPhantoms';
            needSimPhantom = true;
            Mustparam = obj.usParam;

            for k = 1:numFrames
                disp(['--------- Simulation of frame ',num2str(k),' ---------']);
                switch obj.MotionMode
                    case 1 % Only tissue motion
                        load(fullfile(outputFolder,['MovedPhy_',num2str(k),'.mat']));
                        scattersPos = [TissuePoints;obj.flowScatters{k}.pos];
                        scattersRC  = [ones(length(TissuePoints),1);obj.flowScatters{k}.rc];
                    case 2 % Both tissue + vessel shift
                        load(fullfile(outputFolder,['MovedPhy_',num2str(k),'.mat']));
                        deltaD =[2e-3,0,4e-3]./Mustparam.PRF;
                        movedFlow =  obj.flowScatters{k}.pos + deltaD*(k-1);
                        scattersPos = [TissuePoints;movedFlow];
                        scattersRC  = [ones(length(TissuePoints),1);obj.flowScatters{k}.rc];
                    otherwise % No tissue motion
                        scattersPos = [obj.scene.tissueScatters;obj.flowScatters{k}.pos];
                        scattersRC  = [ones(length(obj.scene.tissueScatters),1);obj.flowScatters{k}.rc];
                end

                % Position limits based on probe
                x = scattersPos(:, 1);
                y = scattersPos(:, 2);
                z = scattersPos(:, 3);
                switch obj.globalParam.probeCase
                    case 'L11-5v'
                        x_cond = (x >= -2.1e-2) & (x <= 2.1e-2);
                        y_cond = (y >= -0.5e-2) & (y <= 0.5e-2);
                        z_cond = (z >= 0)   & (z <= 5.1e-2);
                    case 'Vermon'
                        x_cond = (x >= -1.1e-2) & (x <= 1.1e-2);
                        y_cond = (y >= -1.1e-2) & (y <= 1.1e-2);
                        z_cond = (z >= 0)   & (z <= 2.8e-2);
                end
                idx = x_cond & y_cond & z_cond;
                scattersPos = scattersPos(idx, :);
                scattersRC  = scattersRC(idx, :);

                % Define RF output path
                RFoutputName = fullfile('..\Result',obj.globalParam.resultFloder,'RF',['RF_',num2str(k,'%03d'),'.mat']);
                tstart = tic;

                switch obj.MotionMode
                    case {1,2}
                        if obj.globalParam.SimThreeDMode == 1
                            alg.FQSim3D(scattersPos,scattersRC,RFoutputName,struct(obj.usParam));
                        else
                            alg.FQSim2D(scattersPos,scattersRC,RFoutputName,struct(obj.usParam));
                        end
                    otherwise % Combine phantom and vessel simulation
                        if needSimPhantom == true
                            PhantomRF =  alg.FQSim3D_noSave(obj.scene.tissueScatters,...
                                         ones(length(obj.scene.tissueScatters),1),RFoutputName,struct(obj.usParam));
                            needSimPhantom = false;
                        end
                        VesselRF =  alg.FQSim3D_noSave(obj.flowScatters{k}.pos,obj.flowScatters{k}.rc,RFoutputName,struct(obj.usParam));
                        RF = cell(obj.usParam.Na, 1);
                        for kk = 1:length(VesselRF)
                            a = PhantomRF{kk};
                            b = VesselRF{kk};
                            padRows = abs(size(a,1) - size(b,1));
                            if size(a,1) > size(b,1)
                                b_padded = [b; zeros(padRows, size(b, 2), 'like', b)];
                                RF{kk} = a + b_padded;
                            else
                                a_padded = [a; zeros(padRows, size(a, 2), 'like', a)];
                                RF{kk} = b + a_padded;
                            end
                        end
                                   % Save RF data
                folderPath2 = fileparts(RFoutputName);
                if ~isfolder(folderPath2)
                    mkdir(folderPath2);
                end
                save(RFoutputName, 'RF');    
                end

 

                elapsedTime = toc(tstart);
                disp(['Frame time: ', num2str(elapsedTime/60, '%.2f'), ' mins for simulation']);
                disp(['Remaining time (est): ', num2str(elapsedTime*(numFrames-k)/60/60, '%.2f'), ' hours']);
                clear scattersPos scattersRC;
            end
        end
    end
end