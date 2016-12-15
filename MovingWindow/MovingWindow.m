function [dataInWindow, dateWindowEnd, flag] = MovingWindow(stage, data, h, idxMthEnd, iMth)

    idxWindowEnd = idxMthEnd(iMth);
    
    dateWindowEnd = data(idxWindowEnd, h.yyyymmdd);
    
    switch stage
        case 1
            quantityAdded = -2;
            dateField = 'year';
            minObs = 250;
        case 2
            quantityAdded = -1;
            dateField = 'year';
            minObs = 50;
        otherwise
            % pass
    end
    
    
    dateWindowStart = my_addtodate(dateWindowEnd, quantityAdded, dateField);
    
    idxWindowStart = find(data(:, h.yyyymmdd) > dateWindowStart, 1, 'first');
    
    dataInWindow = data(idxWindowStart:idxWindowEnd, :);
    
    if size(dataInWindow, 1) >= minObs
        flag = true;
    else
        flag = false;
    end
      
end