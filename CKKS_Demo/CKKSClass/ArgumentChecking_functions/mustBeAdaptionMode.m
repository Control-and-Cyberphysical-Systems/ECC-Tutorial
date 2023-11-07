function mustBeAdaptionMode(adaption)
    if ~ischar(adaption) 
        error('adaption must be a character vector')
    end
    if ~strcmp(adaption,'auto') &&  ~strcmp(adaption,'manual')
        error("Supported encoding types are 'auto' and 'manual'")
    end
end