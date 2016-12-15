function [fval, Nd1] = ComputeFvalFder(val, data, h, vol, debt, maturity)
    
d = (log(val./debt)+(data(:, h.rfRate)+0.5*vol^2)*maturity)./(vol*sqrt(maturity));

% Nd1 = normcdf(d);
% Nd2 = normcdf(d - vol*sqrt(maturity));
Nd1 = 0.5*erf(d./sqrt(2)) + 0.5;
Nd2 = 0.5*erf((d - vol*sqrt(maturity))./sqrt(2)) + 0.5;
Nd1(Nd1<eps) = eps;
Nd2(Nd2<eps) = eps;
fval = val .* Nd1 - exp(-data(:, h.rfRate)*maturity).*debt.*Nd2;

end