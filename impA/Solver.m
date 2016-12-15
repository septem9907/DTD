function impA = Solver(data, h, vol, debt, maturity, initImpA)

% initialization
if ~isnan(initImpA)
    val = initImpA;
else
    ub = SetUppderBound(data, h, vol, debt, maturity);
    lb = SetLowerBound(data, h, vol, debt, maturity);
    val = 0.5*(ub+lb);
end
% Newton Method
impA = NewtonMethod(val, data, h, vol, debt, maturity);
% Bisection Method(When Newton Method failed)
if isnan(impA)
    ub = SetUppderBound(data, h, vol, debt, maturity);
    lb = SetLowerBound(data, h, vol, debt, maturity);        
    impA = BisectionMethod(data, h, vol, debt, maturity, lb, ub);
end
end