function mustBeKey(key)
    if ~isa(key,'struct')
        error('key must be struct')
    end
    names = fieldnames(key);
    if ~isa(key.(names{1}),'RingElement') || ...
       ~isa(key.(names{2}),'RingElement')
        error('key elements must be RingElements')
    end
end