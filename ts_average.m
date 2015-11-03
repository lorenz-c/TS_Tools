function ts_out = ts_average(ts_in, tres_out, method)
% Compute temporal averages of a time-series.
%--------------------------------------------------------------------------
% Input:        
% - ts_in       Matrix which contains the input time-series. 
% - tres_out    Time scale of the output; can be set to 'annual' (average 
%               over each year), 'ltm' (long term mean), 'mnthly' (long-
%               term monthly mean), 'ssnl' (long-term seasonal mean of 
%               DJF, ...), or 'ssnl_mnthly' (seasonal mean of each year).
% - method      Method for averaging; can be set to 'sum', 'mean', 
%               'sum_squared'; 
% - tres_in     Time scale of ts_in; can be set to 'hourly', 'daily', 
%               'monthly', or 'annual'
% - clms        Ordering of the time-vector; the first, second and third
%               (optionally) element corresponds to the columns containing
%               the years, months, and days, respectively. 
% - miss        Identifier of missing value
% - indx_row    If indx_row == 1, the first row in fld is used as an 
%               indexing row.
%--------------------------------------------------------------------------
% Output:       
% - otpt      	Matrix which contains the desired temporal averages
%
%--------------------------------------------------------------------------
% Author: Christof Lorenz
% Date:   July 2011
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------
% Updates: - 01.04.2015: Brush up code, added some comments
%--------------------------------------------------------------------------

if nargin < 3 | isempty(method), method = 'nanmean'; end
if isfield(ts_in.DataInfo, 'MissingValue')
    miss = ts_in.DataInfo.MissingValue;
else
    miss = NaN;
end

if ~isnan(miss)
    ts_in.Data(ts_in.Data == miss) = NaN;
    miss = NaN;
end


ts_out = ts_in;

ts_out.Data = [];

%--------------------------------------------------------------------------
%                       Annual long term average
%--------------------------------------------------------------------------
if strcmp(tres_out, 'annual_lt') 
    
    ts_out.Data = agg_ts(ts_in.Data, method);

%--------------------------------------------------------------------------
%                      Seasonal long term average
%--------------------------------------------------------------------------
elseif strcmp(tres_out, 'seasonal_lt')
    
    % Create a vector which contains the months where data is available
    mnths = ts_in.DateTime(:, 2);
    
    % Loop over the four seasons
    for i = 1:4
        
        if i == 1
            indx_1 = find(mnths == 12);
            indx_2 = find(mnths == 1);
            indx_3 = find(mnths == 2);            
        else
            indx_1 = find(mnths == i*3-3);
            indx_2 = find(mnths == i*3-2);
            indx_3 = find(mnths == i*3-1);
        end
        
        indx = [indx_1; indx_2; indx_3];
        
        ts_out.Data(i, :) = agg_ts(ts_in.Data(indx, :), method);
        
        clear indx_1 indx_2 indx_3
    end

    ts_out.DateTime = (1:4)';
    ts_out.TimeStamp = (1:4)';
%--------------------------------------------------------------------------
%                     Monthly long term average
%--------------------------------------------------------------------------
elseif strcmp(tres_out, 'monthly_lt')
        
    % Create a vector which contains the months where data is available
    mnths = ts_in.DateTime(:, 2);
    
    % Loop over twelve months 
    for i = 1:12
        % Search for all months in the data which correspond to the ith
        % month:
        indx = find(mnths == i);
        
        ts_out.Data(i, :) = agg_ts(ts_in.Data(indx, :), method);
    end
    
    ts_out.DateTime = (1:12)';
    ts_out.TimeStamp = (1:12)';

%--------------------------------------------------------------------------
%                             Annual average
%--------------------------------------------------------------------------
elseif strcmp(tres_out, 'annual')
    
    % Vector of unique years
    yrs = unique(ts_in.DateTime(:, 1));
    
    for i = 1:length(yrs)
        % Search for all rows which correspond to yrs(i)
        indx = find(ts_in.DateTime(:, 1) == yrs(i));
        ts_out.Data(i, :) = agg_ts(ts_in.Data(indx, :), method);
    end
    
    ts_out.DateTime  = [yrs ones(length(yrs), 2) zeros(length(yrs), 3)];
    ts_out.TimeStamp = datenum(ts_out.DateTime);
%--------------------------------------------------------------------------
%                           Seasonal average
%--------------------------------------------------------------------------
elseif strcmp(tres_out, 'seasonal')
    
    ts_out.DateTime = [];
    ts_out.TimeStamp =[];
    
    % Get the first and last year of the time-series
    syr = ts_in.DateTime(1, 1);
    eyr = ts_in.DateTime(end, 1);
    
    % Vector of unique years
    yrs_unique = unique(ts_in.DateTime(:, 1));
    yrs   = ts_in.DateTime(:, 1);
    mnths = ts_in.DateTime(:, 2);
    
    k = 1;
    for i = 1:length(yrs_unique)
        for j = 1:4
            if j == 1
                indx_1 = find(mnths == 12 & yrs == yrs_unique(i)-1);
                indx_2 = find(mnths == 1  & yrs == yrs_unique(i));
                indx_3 = find(mnths == 2  & yrs == yrs_unique(i));            
            else
                indx_1 = find(mnths == j*3-3 & yrs == yrs_unique(i));
                indx_2 = find(mnths == j*3-2 & yrs == yrs_unique(i));
                indx_3 = find(mnths == j*3-1 & yrs == yrs_unique(i));
            end
            
            indx = [indx_1; indx_2; indx_3];
  
            ts_out.Data(k, :)     = agg_ts(ts_in.Data(indx, :), method);
            ts_out.DateTime(k, :) = [yrs_unique(i) (j-1)*3+1 15 0 0 0];
            
            k = k + 1;
        end
    end
        
    ts_out.TimeStamp = ts_in.TimeStamp(2:3:end, :);

%--------------------------------------------------------------------------
%                           Monthly average
%--------------------------------------------------------------------------
elseif strcmp(tres_out, 'monthly')
    
    ts_out.DateTime = [];
    ts_out.TimeStamp =[];
    
    % Vector of unique years
    yrs = unique(ts_in.DateTime(:, 1));
    
    k = 1;
    % 1. Loop over the years
    for i = 1:length(yrs)
        % 2. Loop over the months
        for j = 1:12
            % Search for all rows which correspond to yrs(i) and mnths j
            indx = find(ts_in.DateTime(:, 1) == yrs(i) & ...
                        ts_in.DateTime(:, 2) == i);
                    
            ts_out.Data(k, :) = agg_ts(ts_in.Data(indx, :), method);
            
            ts_out.DateTime(k, :) = [yrs(i) j 15 0 0 0];
            
            k = k + 1;
        end
    end

    ts_out.TimeStamp = datenum(ts_out.DateTime);
%--------------------------------------------------------------------------
%                          Daily average
%--------------------------------------------------------------------------
elseif strcmp(tres_out, 'daily')
           
    ts_out.DateTime = [];
    ts_out.TimeStamp =[];
    
    % Vector of unique years
    yrs = unique(ts_in.DateTime(:, 1));
    
    k = 1;
    % 1. Loop over the years
    for i = 1:length(yrs) 
        % 2. Loop over the months
        for j = 1:12
            nrd = eomday(yrs(i), j)
            % 3. Loop over each day of a month
            for k = 1:nrd
                % Search for all rows which correspond to yrs(i), mnths j,
                % and days k
                indx = find(ts_in.DateTime(:, 1) == yrs(i) & ...
                            ts_in.DateTime(:, 2) == i     & ...
                            ts_in.DateTime(:, 3) == k);
                    
                ts_out.Data(k, :) = agg_ts(ts_in.Data(indx, :), method);
                
                ts_out.DateTime(k, :) = [yrs(i) j 15 k 0 0];
                
                k = k + 1;
            end
        end   
    end
end


ts_out.DataInfo.TempRes = tres_out;
ts_out.DataInfo.UserInfo = {ts_out.DataInfo.UserInfo; ...
                           [datestr(datetime), ': ', tres_out, ' ', method, ' during ', datestr(ts_in.DateTime(1, :)), ' and ', datestr(ts_in.DateTime(end, :)),]; ...
                           [datestr(datetime), ': ', 'Check the unit of the calculated variable!!!']};


end
            
      
function Data_out = agg_ts(fld, method)
    if strcmp(method, 'mean')          
        Data_out = mean(fld, 1);
    elseif strcmp(method, 'nanmean')
        Data_out = nanmean(fld, 1);
    elseif strcmp(method, 'sum')
        Data_out = sum(fld, 1);
    elseif strcmp(method, 'nansum')
        Data_out = nansum(fld, 1);
    elseif strcmp(method, 'sum_squared')
        Data_out = sum(fld.^2, 1);
    elseif strcmp(method, 'nansum_squared')
        Data_out = nansum(fld.^2, 1); 
    elseif strcmp(method, 'rms')
        Data_out = sqrt(sum(fld.^2, 1));     
    elseif strcmp(method, 'nanrms')
        Data_out = sqrt(nansum(fld.^2, 1)); 
    end
end


















            
        