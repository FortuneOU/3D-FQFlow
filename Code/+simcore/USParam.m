classdef USParam
    properties
        fc       % Center frequency [Hz]
        fs       % Sampling frequency [Hz]
        bandwidth% Fractional bandwidth [%]
        width    % Element width [m]
        height   % Element height [m]
        radius   % Array radius (Inf = linear)
        Nelements% Number of elements
        pitch    % Pitch [m]
        focus    % Focal distance [m]
        PRF      % Pulse repetition frequency [Hz]
        elements % Element positions [2×Nelements]
        fnumber  % f/# for focusing (unused in default)
        c        % Sound speed [m/s]
        t0       % Start time [s]
        passive  % Passive mode flag
        % TX parameters
        Na       % Number of transmissions
        txdel    % Delay laws (3D)
        txdel2d  % Delay laws (2D PW)
    end
    methods
        function obj = USParam(globalParam)
            % Select probe configuration
            switch globalParam.probeCase
                case 'L11-5v'
                    param = [];
                    param.fc = 5.208e6;
                    param.fs = 4*param.fc;
                    param.bandwidth = 100;
                    param.width  = 0.27e-3;
                    param.height = 5e-3;
                    param.radius = Inf;
                    param.Nelements =128;
                    param.pitch = 0.3e-3;
                    param.focus = 0.018;
                    param.PRF   = 800;
                    xe = ((1:param.Nelements)-63.5)*param.pitch;
                    ye = zeros(size(xe));
                    param.elements = [xe(:).'; ye(:).'];

                    % Define TX laws (3D multi-angle)
                    Na = 7;
                    tiltX= -3*pi/180:1*pi/180:3*pi/180;
                    tiltY = zeros(size(tiltX));
                    txdel = cell(Na,1);
                    for k = 1:Na
                        txdel{k} = txdelay3(param,tiltY(k),tiltX(k));
                    end

                    % Define TX for 2D PW
                    tilt = linspace(-3*pi/180,3*pi/180,Na);
                    txdel2d = cell(Na,1);
                    for k = 1:Na
                        txdel2d{k} = txdelay(param,tilt(k));
                    end

                case 'Vermon'
                    load ..\TX\Vermon1024_PWC_SF_RF_sim.mat;
                    param = [];
                    param.fc = 7.81e6;
                    param.bandwidth = 70;
                    param.width  = 270e-6;
                    param.height = 270e-6;
                    param.fs   = 4*param.fc;
                    param.Nelements =1024;
                    param.radius = Inf;
                    param.PRF    = 800;
                    param.pitch  = 0.3e-3;
                    param.fnumber= [0 0];
                    param.c      = 1540;
                    param.t0     = 0;
                    param.passive= false;
                    param.elements = [element_position(:,1)'; element_position(:,2)'];

                    Na = 5;
                    txdel = cell(Na,1);
                    for k = 1:Na
                        txdel{k} = delays(k,:);
                    end
                    txdel2d = zeros(1,1);

                otherwise
                    error('Unknown probeCase: %s',globalParam.probeCase);
            end

            % Map struct fields into object properties
            fields = fieldnames(param);
            for i = 1:numel(fields)
                if isprop(obj, fields{i})
                    obj.(fields{i}) = param.(fields{i});
                end
            end
            obj.Na     = Na;
            obj.txdel  = txdel;
            obj.txdel2d= txdel2d;
        end
    end
end