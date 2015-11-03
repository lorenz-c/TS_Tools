function vars = remdims(inpt)
% The function extracts a list of variables of a datastructure inpt, which
% are not used as dimension variables (such as, e.g., time, lat, lon).
%--------------------------------------------------------------------------
% Input (required):
% - inpt        Datastructure with one or more variables, which are not
%               dimension variables
% Output:
% - vars        List of variables, which are not used as dimension
%               variables
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         November 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------

% Dimensions of the data
dims = inpt.Dimensions;

% Get the names of the variables
vars = fieldnames(inpt.Variables);
    
% Get the variable IDs which do not contain dimension variables
for i = 1:length(dims)
    dim_ids(i) = find(ismember(vars, dims{i}));
end
    
% Remove the dimension variables
vars(dim_ids, :) = [];