function [y,u] = Mod(x,q)
% symmetric modulus [-q/2,q/2)
    arguments
        x
        q {mustBeModulus(q)}
    end
    if isa(x,'RingElement') || isa(x,'Encoding')
        u = floor(x.coefficients/q+vpa(1/2));
        y = x.coefficients-u*q;
    else
        u = floor(x/q+vpa(1/2));
        y = x-u*q;
    end
end