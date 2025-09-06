classdef spaceTransform
    properties
        perm     % Dimension permutation mapping, e.g. [1 2 3], [2 3 1], etc.
        offset1  % First translation offset
        offset2  % Second translation offset
        name     % Optional transform name/identifier
    end
    methods
        % Constructor
        function obj = spaceTransform(perm,offset1,offset2,name)
            if nargin < 4, name = ""; end
            obj.perm    = perm;
            obj.offset1 = offset1;
            obj.offset2 = offset2;
            obj.name    = name;
        end

        % Apply transform to positions (permutation + offsets)
        function out = transform(obj, points)
            temp = points(:, obj.perm);  % dimension reorder
            temp = temp - obj.offset1;   % apply first offset
            temp = temp + obj.offset2;   % apply second offset
            out = temp;
        end

        % Inverse transform (reverse order and apply negative offsets)
        function out = invTransform(obj, points)
            temp = points - obj.offset2;
            temp = temp + obj.offset1;
            % Compute inverse permutation
            invPerm = [0 0 0];
            for i=1:3
                invPerm(obj.perm(i)) = i;
            end
            out = temp(:, invPerm);
        end

        % Apply transform to vectors (only permutation, no offset)
        function out = transformVec(obj, points)
            temp = points(:, obj.perm);
            out = temp;
        end
    end
end