function main()


econID = 15;


addpath([pwd '/DtdGenerator/'])
addpath([pwd '/impA/'])
addpath([pwd '/tools/'])

% load data
try
    data = load([pwd '/../data/input/data_' num2str(econID) '.mat'], 'data');
    data = data.data;
catch
    data = load([pwd '/../data/input/DTDInput_' num2str(econID) '.mat'], 'dtdInput');
    data = data.dtdInput;
    data = PrepareData(data);
    save([pwd '/../data/input/data_' num2str(econID) '.mat'], 'data', '-v7.3')
end

h = struct('cpnyID', 1, 'yyyymmdd', 2, 'mktCap', 3, 'curLiab', 4, ...
    'ltBorrow', 5, 'othLiab', 6, 'totAsset', 7, 'rfRate', 8, 'dateIntvl', 9);

% x = [27731, 31536, 32342, 33438, 34560, 34963, 38545, 39216, 40278, 127029];
% data = data(ismember(floor(data(:, 1)/1000),x), :);

stage = 1;

maturity = 250;


%%
cpnyArray = unique(data(:, 1));
numOfCpny = length(cpnyArray);
output = cell(numOfCpny, 1);


tic
parfor iCpny = 1:numOfCpny
       
    cpnyID = cpnyArray(iCpny);
    
    dtdOutputThisCpny = GenerateDtdOutput4SingleCpny(cpnyID, stage, data, h, maturity);
    
    output{iCpny, 1} = dtdOutputThisCpny;
   
end
toc


%%
output = cell2mat(output);
output(any(~isfinite(output), 2), :) = [];
output = sortrows(output, 1:2);

save([pwd '/../data/output/output_' num2str(econID) '.mat'], 'output')


end