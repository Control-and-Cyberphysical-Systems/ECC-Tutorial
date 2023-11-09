function z = Ecd(x,s,q)
    arguments
        x
        s {mustBePositive}
        q {mustBeModulus(q)}
    end
        z = Mod(round(s*x),q);
end