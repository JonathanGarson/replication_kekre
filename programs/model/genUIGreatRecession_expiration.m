function [dbarPathWeeks,dbarPathMonths,dbarShockMonths] = genUIGreatRecession_expiration(Tbar,lastHoriz)
%% genUIGreatRecession_expiration.m
% Returns path of UI to be used in Great Recession simulation, assuming
%  that benefit extension in 1/2013 was expected for lastHoriz periods
%  but then expires in 12/2013

dbarPathWeeks = 26*ones(Tbar,80);
 % element (i,j) is the perceived level of dbar 'i-1' months 
 % ahead of calendar month 'j', where j=1 corresponds to May 08
dbarShockMonths = zeros(1,80);
 
% July 2008 (calendar month 3)
dbarShockMonths(3) = 1;
dbarPathWeeks(1:12,3) = 26+13;
for time=4:7
	dbarPathWeeks(1:(end-1),time) = dbarPathWeeks(2:end,time-1);
end

% December 2008 (calendar month 8)
dbarShockMonths(8) = 1;
dbarPathWeeks(1:12,8) = 26+33;
for time=9:10
	dbarPathWeeks(1:(end-1),time) = dbarPathWeeks(2:end,time-1);
end

% March 2009 (calendar month 11)
dbarShockMonths(11) = 1;
dbarPathWeeks(1:15,11) = 26+33;
for time=12:18
	dbarPathWeeks(1:(end-1),time) = dbarPathWeeks(2:end,time-1);
end

% May 2009 (calendar month 13) -- assume EB introduced here, given that 
%  Farber-Valletta median across U.S. states exceeds EUC maximum starting now
dbarShockMonths(13) = 1;
dbarPathWeeks(1:13,13) = 26+33+20;
for time=14:19
	dbarPathWeeks(1:(end-1),time) = dbarPathWeeks(2:end,time-1);
end

% November 2009 (calendar month 19)
dbarShockMonths(19) = 1;
dbarPathWeeks(1:7,19) = min(26+53+20,99);
for time=20:20
	dbarPathWeeks(1:(end-1),time) = dbarPathWeeks(2:end,time-1);
end

% January 2010 (calendar month 21)
dbarShockMonths(21) = 1;
dbarPathWeeks(1:7,21) = min(26+53+20,99);
for time=22:22
	dbarPathWeeks(1:(end-1),time) = dbarPathWeeks(2:end,time-1);
end

% March 2010 (calendar month 23)
dbarShockMonths(23) = 1;
dbarPathWeeks(1:6,23) = min(26+53+20,99);

% April 2010 (calendar month 24)
dbarShockMonths(24) = 1;
dbarPathWeeks(1:7,24) = min(26+53+20,99);
for time=25:27
	dbarPathWeeks(1:(end-1),time) = dbarPathWeeks(2:end,time-1);
end

% August 2010 (calendar month 28)
dbarShockMonths(28) = 1;
dbarPathWeeks(1:9,28) = min(26+53+20,99);
for time=29:32
	dbarPathWeeks(1:(end-1),time) = dbarPathWeeks(2:end,time-1);
end

% January 2011 (calendar month 33)
dbarShockMonths(33) = 1;
dbarPathWeeks(1:17,33) = min(26+53+20,99);
for time=34:44
	dbarPathWeeks(1:(end-1),time) = dbarPathWeeks(2:end,time-1);
end

% January 2012 (calendar month 45)
dbarShockMonths(45) = 1;
dbarPathWeeks(1:8,45) = min(26+53+20,99);
for time=46:46
	dbarPathWeeks(1:(end-1),time) = dbarPathWeeks(2:end,time-1);
end

% March 2012 (calendar month 47)
dbarShockMonths(47) = 1;
dbarPathWeeks(1:3,47) = 26+63;
dbarPathWeeks(4:6,47) = 26+53;
dbarPathWeeks(7:10,47) = 26+37; % to match ui_weeks_fv
for time=48:56
	dbarPathWeeks(1:(end-1),time) = dbarPathWeeks(2:end,time-1);
end

% January 2013 (calendar month 57)
dbarShockMonths(57) = 1;
dbarPathWeeks(1:lastHoriz,57) = 26+37; % to match ui_weeks_fv
for time=58:68
	dbarPathWeeks(1:(end-1),time) = dbarPathWeeks(2:end,time-1);
end

% January 2014 (unexpected expiration)
dbarShockMonths(69) = 1;
dbarPathWeeks(:,69) = 26;
for time=70:80
	dbarPathWeeks(1:(end-1),time) = dbarPathWeeks(2:end,time-1);
end
 
dbarPathMonths = round(dbarPathWeeks/4.5);