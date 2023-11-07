function mustBeScalingFactor(d)
    if ~isnumeric(d) && ~isa(d,'sym')
        error('Scaling factor must be numeric')
    end
    if d < 1
        error('Scaling factor must be greater then or equal to 1')
    end
end