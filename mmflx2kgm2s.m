function out = mmflx2kgm2s(inpt)

% The function converts a datastructure with variables in units of
% [mm/month] or [mm/day] to (COARDS and CF-conform) [kg/m2/s]

% First, get the unit of the input data
unit_in = inpt.DataInfo.Unit;

% Check if unit is either mm/month or mm/day
if ~strcmp(unit_in, 'mm/month') & ~strcmp(unit_in, 'mm/day')
    error('Unit of input variable is not [mm/month] or [mm/day]')
end

% Get the time vector
DateTime = inpt.DateTime;

if strcmp(unit_in, 'mm/month')
    % For monthly data, compute the number of days for each month
    nrd = eomday(DateTime(:, 1), DateTime(:, 2));
    
    for i = 1:length(DateTime)
        inpt.Data(i, :) = inpt.Data(i, :)/(nrd(i)*24*60*60);
 
    end
    
elseif strcmp(unit_in, 'mm/day')
    
    for i = 1:length(DateTime)
        for j = 1:size(inpt.Data, 2)
            inpt.Data(i, :) = inpt.Data(i, :)/(24*60*60);
        end
    end
    
end

out                  = inpt;
out.DataInfo.Unit    = 'kg/m^2/s';

if isempty(out.DataInfo.History)
    out.DataInfo.History = [datestr(now), '; mmflx2kgm2s.m: Unit converted from ', unit_in, ' to kg/m2/s.'];
else
    out.DataInfo.History = [out.DataInfo.History; ...
                            [datestr(now), '; mmflx2kgm2s.m: Unit converted from ', unit_in, ' to kg/m2/s.']];
end

    
