function ub = SetUppderBound(data, h, vol, debt, maturity, maxiter)

if nargin == 5
    maxiter = 30;
end

ub = 0.75*(data(:, h.mktCap)+debt);
for i = 1:maxiter
    fval = ComputeFvalFder(ub, data, h, vol, debt, maturity);
    flagValue = fval < data(:, h.mktCap);
    if any(flagValue)
        ub = ub + flagValue.*ub;
    else
        return
    end
end

end