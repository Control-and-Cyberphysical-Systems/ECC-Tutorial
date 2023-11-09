function me = Dec(ct,sk)
% Decrypt using the secret key `sk`
    arguments
        ct
        sk {mustBeA(sk,'RingElement')}
    end
    
    if ~iscell(ct.c0) % "scalar" case
        mustBeA(ct,'Ciphertext');
        sk.q = ct.q;    % ternary (no reduction needed)
        me = ct.c0+ct.c1*sk;
        me.s = ct.s;
        me.type = ct.type;
    else
        for i = 1:size(ct.c0,1)
            for j = 1:size(ct.c0,2)
                mustBeA(ct,'Ciphertext');
                sk.q = ct.q;
                me{i,j} = ct.c0{i,j}+ct.c1{i,j}*sk;
                me{i,j}.s = ct.s;
                me{i,j}.type = ct.type;
            end
        end
    end
end