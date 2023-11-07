function progressbar(c)

    persistent strCR;           %   Carriage return pesistent variable
    strPercentageLength = 10;   %   Length of percentage string (must be >5)
    strDotsMaximum      = 10;   %   The total number of dots in a progress bar

    if isempty(strCR) && ~ischar(c)
        %fprintf('%s','Initialized without string')
        strCR=-1;
    elseif isempty(strCR) && ischar(c)
        fprintf('%s',c);
        strCR = -1;
    elseif ~isempty(strCR) && ischar(c)
        strCR = [];  
        fprintf([c '\n']);
    elseif isnumeric(c)
        c = floor(c);
        percentageOut = [num2str(c) '%%'];
        percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
        nDots = floor(c/100*strDotsMaximum);
        dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
        strOut = [percentageOut dotOut];
        
        if strCR == -1
            fprintf(strOut);
        else
            fprintf([strCR strOut]);
        end
        
        strCR = repmat('\b',1,length(strOut)-1);
    else
        error('Unsupported argument type');
    end
end
