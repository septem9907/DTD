function dtdOutputThisCpny = GenerateDtdOutput4SingleCpny(cpnyID, stage, data, h, maturity)

    
    dataThisCpny = data(data(:, 1) == cpnyID, :);
    
    [mthArray, idxMthEnd] = unique(floor(dataThisCpny(:, 2) / 100), 'last');
    
    initDtd = struct('val', nan, 'vol', nan, 'wgt', nan, 'drift', nan, ...
        'yyyymmdd', nan, 'impA', nan);
    
    dtdOutputThisCpny = nan(length(mthArray), 7, 'double');
    
    for iMth = 1:length(mthArray)
        
        [dataInWindow, dateWindowEnd, flag] = MovingWindow(stage, dataThisCpny, h, idxMthEnd, iMth);
        
        if flag            
            
            [dtd, res] = DtdGenerator(stage, dataInWindow, h, maturity, initDtd);
            
            dtdOutputThisCpny(iMth, :) = [cpnyID, dateWindowEnd, dtd.val, res.exitflag, dtd.drift, dtd.vol, dtd.wgt];
            initDtd = dtd;
        else
            dtdOutputThisCpny(iMth, :) = [cpnyID, dateWindowEnd, nan, nan, nan, nan, nan];
            continue
        end
        
    end
    
end