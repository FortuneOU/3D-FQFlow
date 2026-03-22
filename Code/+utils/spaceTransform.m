classdef spaceTransform
    properties
        perm      % 维度重排, 如 [1,2,3], [2,3,1]等
        offset1   % 第一次平移
        rotation  % [rx, ry, rz] 分别绕 x,y,z 的旋转弧度
        offset2   % 第2次平移
        scale     % 缩放系数 (标量如 2，或 1x3 向量如 [2,2,2])
        name      % 名称可选
    end
    methods
        function obj = spaceTransform(perm, offset1, rotation, offset2, scale, name)
            % 构造函数：增加了 scale 参数
            if nargin < 6, name = ""; end
            if nargin < 5 || isempty(scale), scale = 1; end % 默认缩放系数为1(即不缩放)
            
            obj.perm = perm;
            obj.offset1 = offset1;
            obj.rotation = rotation;
            obj.offset2 = offset2;
            obj.scale = scale;
            obj.name = name;
        end
        
        function R = getRotationMatrix(obj)
            % 生成旋转矩阵 (先绕X, 后绕Y, 最后绕Z)
            rx = obj.rotation(1);
            ry = obj.rotation(2);
            rz = obj.rotation(3);
            Rx = [1 0 0;
                  0 cos(rx) -sin(rx);
                  0 sin(rx) cos(rx)];
            Ry = [cos(ry) 0 sin(ry);
                  0 1 0;
                  -sin(ry) 0 cos(ry)];
            Rz = [cos(rz) -sin(rz) 0;
                  sin(rz)  cos(rz) 0;
                  0 0 1];
            % 总旋转矩阵 (已修复原代码中的 _ 语法错误)
            R = Rz * Ry * Rx;  
        end
        
        function out = transform(obj, points)
            % 一般变换（位置点）
            temp = points(:, obj.perm);     % 维度重排
            temp = temp .* obj.scale;       % 应用缩放系数 (支持标量或1x3向量)
            temp = temp - obj.offset1;      % 第一次平移
            R = obj.getRotationMatrix();    % 旋转
            temp = (R * temp')';            % 对每一行点应用旋转

            temp = temp + obj.offset2;      % 第二次平移
            out = temp;
        end
        
        function out = invTransform(obj, points)
            % 逆变换
            temp = points - obj.offset2;      % 先消除第二次平移
           
            R = obj.getRotationMatrix();
            temp = (R' * temp')';             % 逆旋转
            temp = temp + obj.offset1;        % 恢复第一次平移
             temp = temp ./ obj.scale;         % 消除缩放 (逆缩放)
            % inverse permutation
            invPerm = zeros(1,3);
            for i = 1:3
                invPerm(obj.perm(i)) = i;
            end
            out = temp(:, invPerm);           % 逆维度重排
        end
        
        function out = transformVec(obj, points)
            % 矢量变换（不做平移，只做维度重排、旋转和缩放）
            temp = points(:, obj.perm);
            R = obj.getRotationMatrix();
             temp =  temp  .* obj.scale
            temp = (R * temp')';
            out = temp;          % 应用缩放系数
        end
    end
end
