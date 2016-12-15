function output = my_addtodate(yyyymmdd, quantityAdded, dateField)

    % check yyyymmdd
        % pass

    yyyymmddNum = datenum(num2str(yyyymmdd), 'yyyymmdd');
    output = str2double(datestr(addtodate(yyyymmddNum, quantityAdded, dateField), 'yyyymmdd'));

end