function [x,u] = Dcd(z,s,q)
    arguments
        z
        s {mustBePositive}
        q {mustBeModulus(q)}
    end
    u = round(z/q);
    y = z-u*q;
    x = y/s;
end