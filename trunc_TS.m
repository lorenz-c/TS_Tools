function ts_out = trunc_ts(ts_in, sdte, edte);
% The function truncates a timeseries, which is given as datastructure to
% the period specified by sdte - edte. Note that the function first
% determines the temporal resolution of the input data. Therefore, it is
% assumed that the structure contains a variable "time", which can be
% accessed via ts_in.Data.time. This variable must be a nx6-matrix,
% which contains the following elements (in the exact order!): 
% [yyyy mm dd hh mm ss]. 
% Dates which are not incuded in ts_in are filled with NaNs.
%--------------------------------------------------------------------------
% Input (required):
% - ts_in       CF-conform data structure
% - sdte, edte  Start- and end-date of the desired period. The format can
%               be one of the following options:
%               sdte = yyyy         -> Start on first of Jan. of year yyyy
%               sdte = [yyyy mm]    -> Start on the first day of month mm
%                                      and year yyyy
%               sdte = [yyyy mm dd] -> Start at 00:00 on day dd of month mm
%                                      and year yyyy
%               If sdte = [yyyy mm dd] and the input timeseries contains
%               monthly values, the dd is omitted. 
% Output:
% - out         Output structure which now covers the desired period.
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         November 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: dtevec.m
%--------------------------------------------------------------------------


% Copy the input to the output structure
ts_out = ts_in;

% Temporal Resolution of the input timeseries 
% STILL A WORKAROUND!!!!!!
if ts_in.Data.time(1, 4) ~= ts_in.Data.time(2, 4)     % --> hourly data
    tres = 'hourly';
elseif ts_in.Data.time(1, 3) ~= ts_in.Data.time(2, 3) % --> daily data
    tres = 'daily';
elseif ts_in.Data.time(1, 2) ~= ts_in.Data.time(2, 2) % --> monthly data
    tres = 'monthly';
elseif ts_in.Data.time(1, 1) ~= ts_in.Data.time(2, 1) % --> yearly data
    tres = 'yearly';
end

% Create a vector with all desired time-steps
[DateTime, TimeStamp] = dtevec(sdte, edte, tres);
    
% Length of the output timeseries
ntstps = length(TimeStamp);

% Get a list of variables which are not used as dimension variables
vars = remdims(ts_in);

for i = 1:length(vars)
    % Get the size of the input data
    size_in    = size(ts_in.Data.(vars{i}));
    
    % Set the length of the first dimension (-> time) to match the desired
    % period
    size_in(1) = ntstps;
    
    % Create an empty array for the output time-series
    Data_out = NaN(size_in);

    % Go through each time-step and check if ts_in contains a value
    for j = 1:length(TimeStamp)
        indx = find(TimeStamp(j) == ts_in.TimeStamp);
    
        if ~isempty(indx)
            Data_out(j, :) = ts_in.Data.(vars{i})(indx, :);
        else
            warning(['trunc_ts.m: No values for ', datestr(TimeStamp(j)), '!; Filled with NaN'])
        end
        
        clear indx
    end
    % Write the truncated (or extended) timeseries to the output variable
    ts_out.Data.(vars{i}) = Data_out;

end

% Write the new timesteps to the output variable
ts_out.Data.time = DateTime;
ts_out.TimeStamp = TimeStamp;

% Update the file history
ts_out.DataInfo.history = strcat(ts_out.DataInfo.history, ...
       [datestr(now), '; trunc_ts.m: Selected period from', sdte, ' to ', edte '. \n']);  



