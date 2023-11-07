function mustBeModulus(q)
    if ~isnumeric(q) && ~isa(q,'sym')
        error('Modulus must be numeric')
    end
    if q < 1
        error('Modulus must be greater then or equal to 1')
    end
    if ~(abs(round(q)-q) < eps)
        error('Modulus must be an integer')
    end
end