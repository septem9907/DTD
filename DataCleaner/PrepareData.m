function dataCleaned = PrepareData(data)


    %*****************************************************************
    % input:    n by 8 matrix with the following variables:
    %           1. firm id
    %           2. date
    %           3. market cap
    %           4. current liability
    %           5. long term borrowing
    %           6. total liability
    %           7. total asset
    %           8. risk-free rate
    %           9. interval between 2 valid observations (hn)
    %*****************************************************************

    
    %% remove invalid data
    
    % preallocated a new matrix for output
    dataCleaned = nan( size(data, 1), size(data, 2) + 1, 'double' );
    
    % companyID and trade dates
    dataCleaned(:, 1:2) = data(:, 1:2);
    
    % risk-free rate
    dataCleaned(:, 8) = data(:, 8)/100/250;
    
    % market cap, total liability and total asset should be positive values
    [idx, ~] = find( data(:, 3) > 0 & ...
                     data(:, 6) > 0 & ...
                     data(:, 7) > 0 );
    dataCleaned(idx, [3, 6, 7]) = data(idx, [3, 6, 7]);
    
    % current liability and long-term borrow should not be negative values
    [idx, ~] = find( data(:, 4) >=0 & ...
                     data(:, 5) >=0 );
    dataCleaned(idx, [4, 5]) = data(idx, [4, 5]);
    
    % other liability
    threshold1 = -0.01;
    threshold2 = 0;
    otherLiab = dataCleaned(:, 6) - dataCleaned(:, 4) - dataCleaned(:, 5);    
    flag1 = otherLiab ./ dataCleaned(:, 6) < threshold1;
    dataCleaned(flag1, 6) = nan;
    flag2 = otherLiab < threshold2;
    dataCleaned(flag2&~flag1, 6) = dataCleaned(flag2&~flag1, 4) + dataCleaned(flag2&~flag1, 5);   
    
    dataCleaned(:, 6) = dataCleaned(:, 6) - dataCleaned(:, 4) - dataCleaned(:, 5); 
    
    %% get interval between 2 valid observations of each firm    
    flag = all(isfinite(dataCleaned(:, 3:8)), 2);
    idxValid = find(flag);
    intvl = diff(idxValid);

%     Obs = [idxValid(1); intvl];
%     dataCleaned(idxValid, 9) = Obs;
%     dataCleaned(~flag, :) = [];
%     [~, firstIdx] = unique(dataCleaned(:, 1));
%     dataCleaned(firstIdx, 9) = 0;

    dataCleaned(idxValid, 9) = [1; intvl];
    dataCleaned(~flag, :) = [];
    
    %%
    dataCleaned(dataCleaned<eps) = eps;

      
end