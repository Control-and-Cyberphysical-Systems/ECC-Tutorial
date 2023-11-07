function y = modinv(x,m)
    arguments
        x
        m {mustBeModulus(m)}
    end
    if x < 0
        x = m + mod(x,m);
    end
    t0 = 0;
    t = 1;
    r0 = m;
    r = mod(x,m);
    while r ~= 0
        q = floor(r0/r);
    
        t_tmp = t;
        t = t0-q*t_tmp;
        t0 = t_tmp;
    
        r_tmp = r;
        r = r0-q*r_tmp;
        r0 = r_tmp;
    end
    if r0 > 1
        error('inverse does not exist')
    end
    if t0 < 0
        t0 = t0+m;
    end
    y = t0;
end