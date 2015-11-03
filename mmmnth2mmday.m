function out = mmmnth2mmday(inpt)

% The function converts a datastructure with variables in units of
% [mm/month] to [mm/day] 

% First, get the unit of the input data
unit_in = inpt.DataInfo.Unit;

% Check if unit is either mm/month or mm/day
if ~strcmp(unit_in, 'mm/month') 
    error('Unit of input variable is not [mm/month]')
end

% Compute the number of days for each month
nrd = eomday(inpt.DateTime(:, 1), inpt.DateTime(:, 2));
    
for i = 1:length(nrd)
    for j = 1:size(inpt.Data, 2)
        inpt.Data{i, j} = inpt.Data{i, j}/(nrd(i)*24*60*60);
	end
end
    
out               = inpt;
out.DataInfo.Unit = 'mm/day';

if isempty(out.DataInfo.History)
    out.DataInfo.History = [datestr(now), '; mmmnth2mmday.m: Unit converted from ', unit_in, ' to mm/day.'];
else
    out.DataInfo.History = {out.DataInfo.History; ...
                            [datestr(now), '; mmmnth2mmday.m: Unit converted from ', unit_in, ' to mm/day.']};
end