function strdigits = Decomp(x,q)
    if isa(x,'LWE_Ciphertext')
        lq = log10(x.q);
        c = [x.b, x.A];
    else
        lq = log10(q);
        c = x;
    end
    strdigits = [];
    base = vpa(10);
    for i = 0:lq-1
        Q = c - mod(c, base^(lq-1-i));
        strdigits = [Q/base^(lq-1-i), strdigits];
        c = c - Q;
    end
end