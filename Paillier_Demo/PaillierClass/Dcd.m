function [a,u] = Dcd(x,s,q)
    arguments
        x
        s {mustBePositive}
        q {mustBeModulus(q)}
    end
    u = floor(x/q+vpa(1/2));
    y = x-u*q;
    a = y/s;
end