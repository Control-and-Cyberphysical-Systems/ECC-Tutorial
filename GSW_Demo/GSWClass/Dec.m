% decryption
function me = Dec(ct,sk)
    arguments
        ct
        sk (:,1) {mustBeCoefficients(sk)}
    end

    if isa(ct,'LWE_Ciphertext')
        me = Mod(ct.b+ct.A*sk,ct.q);
    elseif isa(ct,'GSW_Ciphertext')
        if ~iscell(ct.c)
            me = Mod(ct.c(1,:)*[1;sk],ct.q);
        else
            for i = 1:size(ct.c,1)
                for j = 1:size(ct.c,2)
                    me{i,j} = Mod(ct.c{i,j}(1,:)*[1;sk],ct.q);
                end
            end
        end
    end
end