function [combinedStruct] = merge_structs(struct1, struct2)
% merge_structs: Merges the OnlineDecompParameters and Signals of both
% normal and peel off data sets.

% Get field names (assuming both structs have the same fields)
fields = fieldnames(struct1);

% Initialize new struct
combinedStruct = struct();

% Loop through each field and concatenate
for i = 1:numel(fields)
    field = fields{i};
    
    % Get the data from both structs
    data1 = struct1.(field);
    data2 = struct2.(field);
    
    % Concatenate cell arrays
    if iscell(data1) && iscell(data2)
        combinedStruct.(field) = [data1, data2];
    
    % Concatenate row vectors
    elseif isrow(data1) && isrow(data2)
        combinedStruct.(field) = [data1, data2];
    
    % Concatenate column vectors
    elseif iscolumn(data1) && iscolumn(data2)
        combinedStruct.(field) = [data1; data2];
    
    % Concatenate 2D matrices
    elseif ismatrix(data1) && ismatrix(data2)
        if size(data1, 2) == size(data2, 2)  % Same number of columns
            combinedStruct.(field) = [data1; data2];
        elseif size(data1, 1) == size(data2, 1)  % Same number of rows
            combinedStruct.(field) = [data1, data2];
        else
            error(['Matrices in field ' field ' have incompatible dimensions and cannot be concatenated.']);
        end
    
    % Concatenate 3D matrices
    elseif ndims(data1) == 3 && ndims(data2) == 3
        if size(data1, 1) == size(data2, 1) && size(data1, 2) == size(data2, 2)
            combinedStruct.(field) = cat(3, data1, data2);
        elseif size(data1, 2) == size(data2, 2) && size(data1, 3) == size(data2, 3)
            combinedStruct.(field) = cat(1, data1, data2);
        elseif size(data1, 1) == size(data2, 1) && size(data1, 3) == size(data2, 3)
            combinedStruct.(field) = cat(2, data1, data2);
        else
            error(['3D matrices in field ' field ' have incompatible dimensions and cannot be concatenated.']);
        end
    
    % Concatenate vectors (1D arrays)
    elseif isvector(data1) && isvector(data2)
        % Ensure vectors are the same orientation before concatenating
        if (isrow(data1) && isrow(data2)) || (iscolumn(data1) && iscolumn(data2))
            combinedStruct.(field) = [data1, data2];
        else
            error(['Vectors in field ' field ' have different orientations and cannot be concatenated.']);
        end
    
    else
        error(['Unsupported data type for concatenation in field ' field]);
    end
end

end

