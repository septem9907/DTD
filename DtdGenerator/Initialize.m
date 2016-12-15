function initVal = Initialize(stage, data, h, dtd)



if isnan(dtd.vol)
    initVol = 0.02;
else
    initVol = dtd.vol;
%    initVol(initVol<eps) = eps;
end

switch stage
    case 1     
        if isnan(dtd.wgt)
            initWgt = 0.5;
            ubWgt = 1;
            lbWgt = 0;
        else
            initWgt = dtd.wgt;
%            initWgt(initWgt<eps) = eps;
            ubWgt = min((initWgt+0.05),1);
            lbWgt = max((initWgt-0.05),0);
        end
    case 2
        % pass
        initWgt = func();
        ubWgt = [];
        lbWgt = [];
    otherwise
        % pass
end

if isnan(dtd.impA)
    initImpA = nan;
else
    initImpA = nan(size(data, 1), 1, 'double');
        
    [~, lhs, rhs] = intersect(data(:, h.yyyymmdd), dtd.yyyymmdd);
    
    initImpA(lhs) = dtd.impA(rhs);

    initImpA(data(:, h.yyyymmdd) > dtd.yyyymmdd(end)) = dtd.impA(end);
    
    if any(isnan(initImpA))
        initImpA = nan;
    end
    
end


initVal = struct('initVol', initVol, 'initWgt', initWgt, 'ubWgt', ubWgt, 'lbWgt', lbWgt, ...
    'initImpA', initImpA); 


end