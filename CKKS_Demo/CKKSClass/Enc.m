function ct = Enc(setup,m,pk)
    arguments
        setup {mustBeA(setup,'Setup')}
        m
        pk {mustBeKey(pk)}
    end

    if ~iscell(m) % "scalar" case
        if isempty(m.s)
            error('Only encoded values can be encrypted')
        end
        mustBeA(m,'RingElement');
        [c0,c1] = EncScalar(setup,m,pk);
        ct = Ciphertext(setup,m,c0,c1);
    else
        if isempty(m{1,1}.s)
            error('Only encoded values can be encrypted')
        end
        ct = Ciphertext(setup,m{1,1});
        for i = 1:size(m,1)
            for j = 1:size(m,2)
                mustBeA(m{i,j},'RingElement');
                [c0,c1] = EncScalar(setup,m{i,j},pk);
                ct.c0{i,j} = c0;
                ct.c1{i,j} = c1;
            end
        end
    end
end

function [c0,c1] = EncScalar(setup,m,pk)
    r_vec = datasample([-1,0,1],setup.N); % ternary secret
    r = RingElement(vpa(r_vec),setup.Q);
    e0_vec = randn([1,setup.N])*vpa(setup.sigma);
    e0 = RingElement(round(e0_vec),setup.Q);
    e1_vec = randn([1,setup.N])*vpa(setup.sigma);
    e1 = RingElement(round(e1_vec),setup.Q);
    
    m.q = setup.Q;
    M = RingElement(m.coefficients,m.q);
    c0 = r*pk.c0+M+e0;
    c1 = r*pk.c1+e1;
end