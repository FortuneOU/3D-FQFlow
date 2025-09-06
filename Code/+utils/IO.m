
classdef IO
    properties
        data     % Data to read/write
        filename % Target file name
        folder   % Target folder path
    end
    methods
        % Constructor: initialize with filename and optional folder
        function obj = IO(filename, folder)
            if nargin <2
                folder = pwd; % default to current folder
            end
            obj.filename = filename;
            obj.folder   = folder;
            obj.data     = [];
        end

        % Read specific variable from a .mat file
        function obj = readData(obj, varname)
            file_full = fullfile(obj.folder, obj.filename);
            if exist(file_full, 'file')
                S = load(file_full, varname); % selective load
                if isfield(S, varname)
                    obj.data = S.(varname);
                else
                    error("Variable '%s' not found in %s", varname, obj.filename)
                end
            else
                error("File not found: %s", file_full);
            end
        end

        % Write current data to a .mat file under provided variable name
        function writeData(obj, varname)
            if nargin < 2
                varname = 'data'; % use default field name
            end
            S.(varname) = obj.data;
            file_full = fullfile(obj.folder, obj.filename);
            save(file_full, '-struct', 'S'); % save struct fields as variables in mat
        end
    end
end


