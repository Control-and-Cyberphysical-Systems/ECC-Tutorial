function y = Powermod(x,p,N)
    if isa(x,'Ciphertext')
        x = x.c;
    end
    y = 1;
    while p>0
        if mod(p,2) == 1
            y = mod(y.*x,N);
        end
        x = mod(x.*x,N);
        p = floor(p/2);
    end
end