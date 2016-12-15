function dtd = CalEngine(data, h, maturity, vol, wgt, initImpA)

debt = data(:, h.curLiab)+0.5*data(:, h.ltBorrow)+wgt*data(:, h.othLiab);

impA = Solver(data, h, vol, debt, maturity, initImpA);

logVA = log(impA./data(:, h.totAsset));
logVAdif = diff(logVA);
drift = sum(logVAdif)/sum(data(2:end, h.dateIntvl));

val = (log(impA(end)/debt(end))+drift*maturity)/(vol*sqrt(maturity));

dtd = struct('val', val, 'vol', vol, 'wgt', wgt, 'drift', drift, ...
    'yyyymmdd', data(:, h.yyyymmdd), 'impA', impA);

end