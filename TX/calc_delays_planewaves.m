function [delays_dis, normal_vector] = calc_delays_planewaves(alpha, theta, element_position)
    % 计算多角度平面波的延时
    % 参数:
    %   alpha: 沿 x 轴的旋转角度
    %   theta: 沿 y 轴的旋转角度
    %   element_position: 阵元的位置 (no_element x 3)
    % 返回:
    %   delays: 各个角度下的阵元延时 (nbeams x no_element)
    %   normal_vector: 法向量 (nbeams x 3)

    no_element = size(element_position, 1);  % 阵元数
    nbeams = length(alpha);  % 平面波角度数

    delays_dis = zeros(nbeams, no_element);  % 延时矩阵
    normal_vector = zeros(nbeams, 3);  % 法向量矩阵

    for index = 1:nbeams
        azimuth = alpha(index);
        elevation = theta(index);

        % 将角度转换为弧度
        azimuth = deg2rad(90 - azimuth);
        elevation = deg2rad(90 - elevation);

        % 计算单位向量的分量
        x_component = cos(azimuth);
        y_component = cos(elevation) * sin(azimuth);
        z_component = sin(elevation) * sin(azimuth);

        % 创建并归一化向量
        vector = [x_component, y_component, z_component];
        vector = vector / norm(vector);

        % 存储法向量
        normal_vector(index, :) = vector;

        % 计算延时，并减去最小值进行归一化
        delays_dis(index, :) = element_position * vector' - min(element_position * vector');
    end
end