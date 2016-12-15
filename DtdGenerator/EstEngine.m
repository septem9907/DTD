function res = EstEngine(stage, data, h, maturity, initVal, adjVol, adjWgt, tol)


if nargin == 5
    adjVol = 0.01;
    adjWgt = 1;
    tol = 1e-9;
end

%% prepare input for estimation

switch stage
    case 1
        x0 = [initVal.initVol/adjVol; initVal.initWgt/adjWgt];
        ub = [Inf, initVal.ubWgt/adjWgt];
        lb = [0, initVal.lbWgt/adjWgt];
    case 2
        x0 = initVal.initVol/adjVol;
        ub = Inf;
        lb = 0;
    otherwise
        % pass
end
        
%% estimation process

% trust-region-reflective algorithm
options = optimset('Algorithm','trust-region-reflective','display','notify','TolFun',tol,'TolCon',tol,'GradObj','on','Hessian','on');
[x, fval, exitflag, output] = fmincon(@(x) LikelihoodFunc(x,stage,data,h,maturity,initVal,adjVol,adjWgt),x0,[],[],[],[],lb,ub,[],options);

if output.firstorderopt > 1e-4
    % sqp algorithm 
    options = optimset('Algorithm','sqp','display','notify','TolFun',tol,'TolCon',tol,'GradObj','on');
    [x_cand, fval_cand, exitflag_cand, output_cand] = fmincon(@(x) LikelihoodFunc(x,stage,data,h,maturity,initVal,adjVol,adjWgt),x0,[],[],[],[],lb,ub,[],options);
    
    if fval_cand < fval
        x = x_cand;
        fval = fval_cand;
        exitflag = exitflag_cand;
    end
    
end


%% result
switch stage
    case 1
        res = struct('x', [adjVol*x(1); adjWgt*x(2)], 'fval', fval, 'exitflag', exitflag);
    case 2
        res = struct('x', [adjVol*x; initVal.initWgt], 'fval', fval, 'exitflag', exitflag);
    otherwise
        % pass
end

end