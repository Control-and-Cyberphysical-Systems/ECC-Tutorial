function [y,u] = Mod(z,q)
    arguments
        z
        q {mustBeModulus(q)}
    end
    if isa(z,'RingElement') || isa(z,'Encoding')
        u = round(z.coefficients/q);
        y = z.coefficients-u*q;
    else
        u = round(z/q);
        y = z-u*q;
    end
end