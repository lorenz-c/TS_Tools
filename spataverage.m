function otpt = spataverage(inpt, id_map, area_id, varargin)
% The function computes time-series of area-weighted means over selected 
% areas. These areas are defined by an id_map (a matrix where connected
% regions have the same id) and an area_id vector (or scalar), which 
% contains all the areas over which the input fields should be
% aggregated
%--------------------------------------------------------------------------
% Input:    - inpt      Cell array or structure which contains the input 
%                       fields.                
%           - id_map    Map which defines the different areas
%           - area_id   Vector (or scalar) which contains the ids of the 
%                       desired areas
%
%           - theta     Latitudes of the gridcells in inpt 
%                       default: 89.75°:-0.5°:-89.75°;
%           - dlambda   Size of a gridcell in longitudinal direction
%                       default: 0.5°
%           - clms      Columns which contain the time data 
%                       (in y,m,d,h,...). For cell-arrays: the last element
%                       of clms must point on the column which contains the
%                       actual data
%           - miss      Identifier for missing values in the input fields
%                       default: NaN
%           - method    The function computes either a weighted mean
%                       (wmean), regular mean (mean), weighted sum (sum), 
%                       or regular sum (sum) over the regions defined in 
%                       the id_map
%                       default: wmean
%           - gridarea  Method for calculating the size of the gridcells
%                       Can be set to regular (default), haversine, cos, or
%                       vincenty
%           - tres_in   Defines the temporal resolution of the input
%                       dataset. Can be set to hourly, daily, monthly 
%                       (default), or annual
%           - squared   If set to true, the output will be squared 
%                       default: false
%           - frmt_out  Defines the output format. Can be set to matrix
%                       or struct (default)
%                                      
% Output:   - otpt      Matrix or structure array, which contains the 
%                       aggregated values over the regions which were 
%                       defined by the id_map and the area_ids.
%--------------------------------------------------------------------------
% Author: Christof Lorenz
% Date:   July 2011
%--------------------------------------------------------------------------
% Uses: area_wghts.m, cell2catchmat.m, copy_datastruct.m,
% create_datastruct.m, DataConvert.m
%--------------------------------------------------------------------------
% Updates: - 25.01.2013 For-loops removed, switched to matrix-based 
%                       computation 
%          - 14.10.2015 Added support for new data structures, brushed up
%                       help text and code
%--------------------------------------------------------------------------

% Checking input arguments and setting default values
pp = inputParser;
pp.addRequired('inpt', @(x) (iscell(x) | isstruct(x) | ismatrix(x)));   
pp.addRequired('id_map', @(x) (isnumeric(x) | iscell(x)) | isstruct(x));
pp.addRequired('area_id', @isnumeric);                     

pp.addParamValue('theta', (89.75:-0.5:-89.75)', @isnumeric);
pp.addParamValue('dlambda', 0.5, @isnumeric);
pp.addParamValue('clms', [ ], @isnumeric);              
pp.addParamValue('miss', -9999, @(x) (isnumeric(x) | strcmp(x, 'NaN')));
pp.addParamValue('method', 'wmean')
pp.addParamValue('gridarea', 'regular')
pp.addParamValue('tres_in', 'monthly')
pp.addParamValue('squared', false)
pp.addParamValue('frmt_out', 'struct')

pp.parse(inpt, id_map, area_id, varargin{:})

clms       = pp.Results.clms;
miss       = pp.Results.miss;
method     = pp.Results.method;
theta      = pp.Results.theta;
dlambda    = pp.Results.dlambda;
gridarea   = pp.Results.gridarea;
tres_in    = pp.Results.tres_in;
squared    = pp.Results.squared;
frmt_out   = pp.Results.frmt_out;

clear pp

if isstruct(inpt)
    dataformat = 1;
elseif iscell(inpt)
    dataformat = 2;
elseif ismatrix(inpt) | ~iscell(inpt)
    dataformat = 3;
end

% First, separate data and time information from the input data
[Data, DateVector] = DataConvert(inpt, tres_in, frmt_out, clms);

% Check if the input data structure contains a field for missing values
if dataformat == 1
    if isfield(inpt, 'Info')
        if isfield(inpt.Info, 'mval')
            miss = inpt.Info.mval;
        end
    elseif isfield(inpt, 'info')
        if isfield(inpt.info, 'mval')
            miss = inpt.info.mval;
        end
    end
end

% Number of regions and time-steps
nr_catch = length(area_id);  
nts      = size(Data, 1);

if isstruct(id_map)
    mask_info = id_map.MaskInfo;
    id_map    = id_map.Data;
else
    mask_info = [];
end

% Create a binary mask to remove all "unwanted" elements
bin_mask = zeros(size(Data{1}));
for i = 1:nr_catch
    bin_mask(id_map == area_id(i)) = 1; 
end



% Go through the input dataset and set all missing elements to zero in the
% binary and the id map
for i = 1:nts
    if isnan(miss)
        bin_mask(isnan(Data{i}))  = 0;
    else
        bin_mask(Data{i} == miss) = 0;
    end
end

id_map(bin_mask == 0) = 0;


% Re-arrange the maps to vectors which contain only the non-zero elements
% of bin_mask
bin_vec = bin_mask(bin_mask ~= 0);
id_vec  = id_map(bin_mask ~= 0);

% For each region, a row in the matrix H is created. At this stage, H is
% binary and defines if an element in the data matrix contains to the
% current catchment or not.
for i = 1:nr_catch
    tmp = bin_vec;
    tmp(id_vec ~= area_id(i)) = 0;
    H(:, i) = tmp;
end

% For weighted computations, a matrix A_mer is created which contains the
% areas of the pixels. This is used, depending on the chosen method, to
% apply area weights to the elements in the data matrix. As these are all 
% linear operations (y=A*H), the weights are added to the H-matrix.
if strcmp(method, 'wmean')   
    
    A_mer = area_wghts(theta', dlambda, 'mat', gridarea);
    A_mer = A_mer(bin_mask ~= 0);

    H = H.*repmat(A_mer, [1 nr_catch]);
    H = H./(ones(length(bin_vec), 1)*sum(H));
    
elseif strcmp(method, 'mean')   
    
    H = H./(ones(length(bin_vec), 1)*sum(H));
    
elseif strcmp(method, 'wsum')    
    A_mer = area_wghts(theta', dlambda, 'mat', gridarea);
    A_mer = A_mer(bin_mask ~= 0);
    
    H = H.*repmat(A_mer, [1 nr_catch]);
    
else
    if ~strcmp(method, 'sum')
        error('Method not defined!');
    end
    
end

% Re-arrange the input fields in a big matrix, which contains only the
% pixels which are located in one of the areas of interest
Data_mat = cell2catchmat(Data, bin_mask);

% Now, the aggregated values can be computed by simple matrix
% multiplication:
if squared == true
    Out_Data = sqrt((Data_mat.^2)*(H.^2));
else
    Out_Data = Data_mat*H;
end
    
% Create the output either as a matrix or as a structure
if strcmp(frmt_out, 'struct')
    if dataformat == 1
        otpt                       = copy_datastruct(inpt, 'spatial', 'aggregated');
    elseif dataformat == 2
        otpt                       = create_datastruct('aggregated');
        otpt.DataInfo.Center       = inpt{1, 1};
        otpt.DataInfo.Unit         = inpt{1, end};
        otpt.DataInfo.VariableName = inpt{1, 2};
    elseif dataformat == 3
        otpt = create_datastruct('aggregated');
    end
        
    otpt.DateTime  = DateVector.DateTime;
    otpt.TimeStamp = DateVector.TimeStamp;
    
    otpt.Aggregate.Mask = id_map;
    otpt.Aggregate.IDs  = area_id;
    if ~isempty(mask_info)
        otpt.Aggregate.IDNames = mask_info.IDNames;
    end
    
    otpt.Data           = Out_Data;
    
    otpt.DataInfo.Type  = 'aggregated';
    otpt.DataInfo.Owner = 'C. Lorenz';

    if isfield(otpt.DataInfo, 'History')
       otpt.DataInfo.History = {otpt.DataInfo.History; ...
                                   [datestr(datetime), ' spataverage.m: computed spatial', method]};
    else
        otpt.DataInfo.History = [datestr(datetime), ' spataverage.m: computed spatial ', method];
    end
      
    
else
    otpt = [zeros(1, 6)         0                    area_id; ...
            DateVector.DateTime DateVector.TimeStamp Out_Data];
end

    
    
    
        
        
        
       
    
    
        
    


