%% estimateSepRate.m
% Estimates separation rate adjusted for time aggregation

tic;

%% Read in dataset

df = readtable('Output/seprate_for_matlab.csv');
df(df.year < 1948, :) = [];
df(df.year > 2019, :) = [];

year = df.year;
month = df.month;
unemploy = df.unemploy;
employ = df.employ;
short_ue_adj = df.short_ue_adj;

%% Compute separation rate without adjusting for time aggregation

seprate = short_ue_adj(2:end)./employ(1:end-1);

%% Compute separation rate adjusting for time aggregation

F = 1 - (unemploy(2:end)-short_ue_adj(2:end))./unemploy(1:end-1);
f = -log(1-F);
nonlineq = @(x) unemploy(2:end) - ...
    ((1-exp(-f-x)).*x./(f+x)).*(employ(1:end-1)+unemploy(1:end-1)) - ...
    exp(-f-x).*unemploy(1:end-1);
seprate_timeagg = fsolve(nonlineq,seprate);

%% Fit linear trend line

start1990 = find(df.year == 1990 & df.month == 1);
cov_1 = cov([[1:(length(seprate_timeagg)-start1990+1)]',seprate_timeagg(start1990:end)]);
slope = cov_1(1,2)/cov_1(1,1);
intercept = mean(seprate_timeagg(start1990:end) - slope*[1:(length(seprate_timeagg)-start1990+1)]');
bestfit = intercept + slope*[1:(length(seprate_timeagg)-start1990+1)]';

%% Characterize detrended separation rate over 1990-2019 and save to file

seprate_timeagg_dt = seprate_timeagg(start1990:end) - bestfit;

cov_2 = cov([seprate_timeagg_dt(1:end-1),seprate_timeagg_dt(2:end)]);
ar1 = cov_2(1,2)/cov_2(1,1);
seprate_timeagg_dt_inno = seprate_timeagg_dt(2:end) - ar1*seprate_timeagg_dt(1:end-1);
disp(['AR(1) detrended, 90-19 = ',num2str(ar1)]);
csvwrite('Output/seprate.csv',...
    [year(start1990:end-1),month(start1990:end-1),seprate_timeagg(start1990:end),seprate_timeagg_dt,[0;seprate_timeagg_dt_inno]]);
% Columns: year, month, separation rate (level), separation rate (detrended), separation rate (innovation)

disp(['Code completed in ',num2str(toc),' s.']);