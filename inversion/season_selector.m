function time_range = season_selector(filename, season, start_year, end_year)
% returns the indices of time instances that belong to a given season

 error(nargchk(4, 4, nargin))

 time  = timenc(filename);

 month = time(:,2);
 year  = time(:,1);

 min_time = min(find(year>=start_year));
 max_time = max(find(year<=end_year));

 month_range = month(min_time:max_time);

 if findstr('ann', season)
  time_range = [1:length(month_range)];
 elseif findstr('djf', season)
  time_range = find(month_range==12|month_range==1|month_range==2)';
 elseif findstr('mam', season)
  time_range = find(month_range==3|month_range==4|month_range==5)';
 elseif findstr('jja', season)
  time_range = find(month_range==6|month_range==7|month_range==8)';
 elseif findstr('son', season)
  time_range = find(month_range==9|month_range==10|month_range==11)';
 elseif findstr('jan', season)
  time_range = find(month_range==1)';
 elseif findstr('feb', season)
  time_range = find(month_range==2)';
 elseif findstr('mar', season)
  time_range = find(month_range==3)';
 elseif findstr('apr', season)
  time_range = find(month_range==4)';
 elseif findstr('may', season)
  time_range = find(month_range==5)';
 elseif findstr('jun', season)
  time_range = find(month_range==6)';
 elseif findstr('jul', season)
  time_range = find(month_range==7)';
 elseif findstr('aug', season)
  time_range = find(month_range==8)';
 elseif findstr('sep', season)
  time_range = find(month_range==9)';
 elseif findstr('oct', season)
  time_range = find(month_range==10)';
 elseif findstr('nov', season)
  time_range = find(month_range==11)';
 elseif findstr('dec', season)
  time_range = find(month_range==12)';
 else
  error(['season not found: ', season])
 end


 time_range = time_range+min_time-1;



