function [dtd, res] = DtdGenerator(stage, data, h, maturity, dtd)


%% initialization
initVal = Initialize(stage, data, h, dtd);


%% estimation

res = EstEngine(stage, data, h, maturity, initVal);


%%
vol = res.x(1);
wgt = res.x(2);
initImpA = initVal.initImpA;


%% calculation

dtd = CalEngine(data, h, maturity, vol, wgt, initImpA);


end