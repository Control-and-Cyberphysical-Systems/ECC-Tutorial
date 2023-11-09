function m = Dec(ct,sk)
    % decryption
    arguments
        ct {mustBeA(ct,'Ciphertext')}
        sk
    end
    skinv = modinv(sk,ct.pk);
    tmp = (Powermod(ct.c,sk,ct.q)-1)./ct.pk;
    m = mod(tmp.*skinv,ct.pk);
end