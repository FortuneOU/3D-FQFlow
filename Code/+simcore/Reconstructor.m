classdef Reconstructor
    properties
        ReconArea     % Imaging reconstruction grid
        ReconRFPath   % Folder containing recorded RF data
        usParam       % Ultrasound probe parameters
        globalParam   % Global parameter set
    end
    methods
        % Constructor
        function obj = Reconstructor(globalParam,usParam)
            obj.globalParam = globalParam;
            obj.usParam     = usParam;
        end

        % Main reconstruction step (2D or 3D beamforming)
        function [] = runRecon(obj)
            floderPath = fullfile('..\Result',obj.globalParam.resultFloder);
            files = dir(fullfile([floderPath,'/RF'],'**','*'));        
            files = files(~[files.isdir]);                     
            numRF = numel(files);

            if obj.globalParam.ReconThreeDMode == 0
                % --- 2D beamforming reconstruction ---
                for k = 1:numRF
                    RFoutputName = fullfile(floderPath,'RF',['RF_',num2str(k,'%03d'),'.mat']);
                    IQoutputName = fullfile(floderPath,'IQ',['IQ_',num2str(k,'%03d'),'.mat']);
                    tstart = tic;
                    alg.FQBmFr2D(RFoutputName,IQoutputName,struct(obj.usParam),obj.globalParam.ReconArea);
                    elapsedTime = toc(tstart);
                    disp(['Frame completed in ', num2str(elapsedTime/60, '%.2f'), ' mins']);
                    disp(['Remaining time (est): ', num2str(elapsedTime*(numRF-k)/60/60, '%.2f'), ' hours']);
                end
            else
                % --- 3D beamforming reconstruction ---
                RFoutputFloder = fullfile(floderPath,'RF');
                IQoutputFloder = fullfile(floderPath,'IQ');
                tstart = tic;
                alg.FQBmFr3D(RFoutputFloder,IQoutputFloder,struct(obj.usParam),obj.globalParam.ReconArea);
                elapsedTime = toc(tstart);
                disp(['3D DAS completed in ', num2str(elapsedTime/60, '%.2f'), ' mins']);
            end
        end
    end
end