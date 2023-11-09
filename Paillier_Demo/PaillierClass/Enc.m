function ct = Enc(m,pk)
    arguments
        m
        pk {mustBeModulus(pk)}
    end
    q = pk^2;
    ct = Ciphertext(vpa([]),pk,q);
    for i = 1:size(m,1) 
        for j = 1:size(m,2) 
            tmp = EncScalar(m(i,j),pk,q);
            ct.c(i,j) = tmp.c;
        end
    end
end


function ct = EncScalar(m,pk,q)
    % encryption routine for a scalar input
    r = vpa(floor(rand()*pk));
    while gcd(r,pk-1) ~= 1
        r = vpa(floor(rand()*pk)); % r\in (Z/pkZ)*; not Pseudorandom!
    end
    tmp1 = Powermod(pk+1,m,q);
    tmp2 = Powermod(r,pk,q);       % slower than matlab's but necessary
    ct = mod(tmp1*tmp2,q);
    ct = Ciphertext(ct,pk,q);
end