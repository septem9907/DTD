function val = NewtonMethod(val, data, h, vol, debt, maturity, tol, maxiter)

if nargin == 6
    tol = 1e-8;
    maxiter = 10;
elseif nargin == 7
    maxiter = 10;
end

for i = 1:maxiter
    [fval, fder] = ComputeFvalFder(val, data, h, vol, debt, maturity);
    step = -(fval - data(:, h.mktCap))./fder;
    if max(abs(step./val)) < tol
        return
    else
        val = val+step;
    end
end

val = nan;

end