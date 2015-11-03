function out = mmflx2kgm2s(inpt, type, varargin)
% The function converts a datastructure with flux variables in units of
% [kg/m2/s] to [mm/month] or [mm/day]. Please note that for the conversion
% to [mm/month], the function requires the respective month of each
% time-step. Therefore, it is assumed that the input structure contains a
% variable "time", where the dates are stored in 6 columns 
% (YYYY MM DD HH MM SS), from which the first two contain the years and
% months.
%--------------------------------------------------------------------------
% Input (required):
% - inpt        CF-conform data structure
% - type        Defines the unit of the output. Possible selections are
%                   1: Convert to [mm/month]
%                   2: Convert to [mm/day] (default)
% - varargin    Names of the variables which should be converted; if
%               varargin is empty, the function transforms all variables
%               (except the dimension variables time, lat, lon, and some
%               z-axis) to the new units
%
% Output:
% - out     
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         November 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------

% If type is empty, set the output unit to [mm/day]
if isempty(type)
    type = 2;
end

% Simply copy the input to the output variable
out = inpt;

if isempty(varargin)
    % Get a list of variables which are not used as dimension variables
    vars = remdims(ts_in);
else
    vars = varargin;
end

% Check the unit of the variables
if type == 1
    % Conversion to [mm/month]
    for i = 1:length(vars)
        if ~strcmp(inpt.Variables.(vars{i}).units, 'kg/m^2/s')
            warning(['mmflx2kgm2s.m: Variable ', vars{i}, ' is not in units of [kg/m^2/s]; skipping variable!']);
        else
            
            % Get the years and months from the input data
            DateTime = inpt.Data.time;
            Yrs      = DateTime(:, 1);
            Mnths    = DateTime(:, 2);
            
            % Compute the number of days for each month
            nrd      = eomday(Yrs, Mnths);
            
            % Get the size of the data
            repsize    = size(inpt.Data.(vars{i}));
            repsize(1) = 1;
            
            % Re-shape the nrd-vector to mach the size of the data
            nrd      = repmat(nrd, repsize)*24*3600;
            
            % Multiply each data slice with the corresponding number of
            % days (x 60 seconds x 60 minutes x 24 hours)
            out.Data.(vars{i}) = inpt.Data.(vars{i}).*nrd;
            
            % Update the MetaData
            out.Variables.(vars{i}).units     = 'mm/month';
            out.Variables.(vars{i}).valid_max = ...
                                 max(out.Data.(vars{i})(:), [], 'omitnan');
            out.Variables.(vars{i}).valid_min = ...
                                 min(out.Data.(vars{i})(:), [], 'omitnan');
        end
    end
    
    % Update the file history
    out.DataInfo.history = strcat(out.DataInfo.history, ...
                    [datestr(now), '; kgm2s2mmflx.m: Converted units from [kg/m^2/s] to [mm/month]. \n']);    

elseif type == 2
    
    % Conversion to [mm/day]
    for i = 1:length(vars)
        if ~strcmp(inpt.Variables.(vars{i}).units, 'kg/m^2/s')
            warning(['mmflx2kgm2s.m: Variable ', vars{i}, 'is not in units of [kg/m^2/s]; skipping variable!']);
        else
            
            % Multiply each data slice with 60 seconds x 60 minutes x 24
            % hours
            out.Data.(vars{i})                = inpt.Data.(vars{i})*3600*24;
            
            % Update the MetaData
            out.Variables.(vars{i}).units     = 'mm/day';
            out.Variables.(vars{i}).valid_max = ...
                                 max(out.Data.(vars{i})(:), [], 'omitnan');
            out.Variables.(vars{i}).valid_min = ...
                                 min(out.Data.(vars{i})(:), [], 'omitnan');
        end
    end
    
    % Update the file history   
    out.DataInfo.history = strcat(out.DataInfo.history, ...
        [datestr(now), '; kgm2s2mmflx.m: Converted units from [kg/m^2/s] to [mm/day]. \n']);    
end
           
