function mustBeEncoder(Encoder)
    if ~isempty(Encoder) && ~isa(Encoder,'CoefficientEncoder')
        error('Encoder must be empty or an CoefficientEncoder')
    end
end