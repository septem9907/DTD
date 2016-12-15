function val = BisectionMethod(data, h, vol, debt, maturity, lb, ub, tol)
if nargin == 7
    tol = 1e-8;
end
maxiter = ceil(log(tol)/log(0.5));
for i = 1:maxiter
    val = 0.5*(ub+lb);
    fval = ComputeFvalFder(val, data, h, vol, debt, maturity);
    flag = fval > data(:, h.mktCap);
    ub = flag.*val + (1-flag).*ub;
    lb = flag.*lb + (1-flag).*val;
    if max((ub-lb)/lb) <= tol
        val = 0.5*(ub+lb);
        return
    end
    val = 0.5*(ub+lb);
end
end