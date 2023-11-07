function mustBeEncodingType(type)
    if isempty(type)
        return
    end
    if ~ischar(type) 
        error('Type must be a character vector')
    end
    if ~strcmp(type,'scalar') &&  ~strcmp(type,'packed')
        error("Supported encoding types are 'scalar' and 'packed'")
    end
end