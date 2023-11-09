function mustBeCtTuple(c)
    if ~iscell(c)
        if ~isempty(c) && ~isa(c,'RingElement')
            error('ciphertext tuple must be empty or RingElement')
        end
    else
        for i = 1:size(c,1)
            for j = 1:size(c,2)
                if ~isempty(c{i,j}) && ~isa(c{i,j},'RingElement')
                    error('ciphertext tuple must be empty or RingElement')
                end
            end
        end
    end
end