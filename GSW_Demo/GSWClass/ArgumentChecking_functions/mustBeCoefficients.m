function mustBeCoefficients(coefficients)
    try
        coefficients = coefficients{1,1};
    end
    if ~isnumeric(coefficients) && ~isa(coefficients,'sym')
        error('Coefficients must be numeric')
    end
end