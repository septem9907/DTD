function lb = SetLowerBound(data, h, vol, debt, maturity, maxiter)

if nargin == 5
    maxiter = 30;
end

lb = 0.6*(data(:, h.mktCap)+debt);
for i = 1:maxiter
    fval = ComputeFvalFder(lb, data, h, vol, debt, maturity);
    flagValue = fval > data(:, h.mktCap);
    if any(flagValue)
        lb = lb - 0.5*flagValue.*lb;
    else
        return
    end
end

end