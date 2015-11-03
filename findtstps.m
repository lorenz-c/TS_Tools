function ts_out = trunc_ts(ts_in, sdte, edte);


% Create a vector with all desired time-steps
[DateTime, TimeStamp] = dtevec(sdte, edte, tres);



% Temporal Resolution of the input timeseries
tres = ts_in.DataInfo.TempRes;



% Length of the output timeseries
ntstps = length(TimeStamp);

% Dimensions of the data
dims = size(ts_in.Data);
    
% Create an empty array for the output time-series
Data_out = NaN([ntstps dims(2:end)]);

% Go through each time-step and check if ts_in contains a value
for i = 1:length(TimeStamp)
    indx = find(TimeStamp(i) == ts_in.TimeStamp);
    
    if ~isempty(indx)
        Data_out(i, :) = ts_in.Data(indx, :);
    else
        warning(['No values for ', datestr(TimeStamp(i)), '!; Filled with NaN'])
    end
end

ts_out           = ts_in;
ts_out.Data      = Data_out;
ts_out.DateTime  = DateTime;
ts_out.TimeStamp = TimeStamp;

ts_out.DataInfo.History = [ts_out.DataInfo.History; [datestr(now), '; findtstps.m: changed period to ', datestr(TimeStamp(1)), ' - ', datestr(TimeStamp(end))]];



