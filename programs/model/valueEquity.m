function q = valueEquity(m,div,qTplus1)
%% valueEquity.m
% Returns (end of period) value of equity q at each date 1,...,T given
%   inverse real interest rates m between dates 1,...,T and dates 2,...,T+1
%   dividends div at dates 1,...T+1
%   terminal value qTplus1 at date T+1

T = length(m);
q = nan(1,T+1);
q(T+1) = qTplus1;
for t=T:-1:1
    q(t) = m(t)*(div(t+1)+q(t+1));
end

q = q(1:T);
