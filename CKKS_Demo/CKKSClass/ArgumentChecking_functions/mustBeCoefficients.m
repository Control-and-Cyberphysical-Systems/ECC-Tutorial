function mustBeCoefficients(coefficients)
    if ~isnumeric(coefficients) && ~isa(coefficients,'sym')
        error('Coefficients must be numeric')
    end
end