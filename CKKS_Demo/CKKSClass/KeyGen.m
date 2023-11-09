function [pk,sk,setup] = KeyGen(setup)
    % Generate the public, secret and evaluation key. The latter is
    % stored as a property of the returned `CKKS` instance `setup`.
    arguments
        setup {mustBeA(setup,'Setup')}
    end
    % NOTE: The randomness and emulated distributions do not satisfy
    %       cryptographic requirements (see comments below)!
    
    % secret
    s_vec = datasample([-1,0,1],setup.N);  % ternary secret
    secret = RingElement(vpa(s_vec),setup.Q);
    sk = secret;
    
    % public key
    a_vec = vpa(rand([1,setup.N]))*setup.Q-setup.Q/2; % uniform
    a = RingElement(round(a_vec),setup.Q);
    e_vec = randn([1,setup.N])*vpa(setup.sigma);      % Gaussian
    e = RingElement(round(e_vec),setup.Q);
    b = -a*secret+e;
    
    pk.c0 = b;
    pk.c1 = a;
    
    % evaluation key
    PQ = setup.P*setup.Q;
    a_vec = vpa(rand([1,setup.N]))*PQ-PQ/2;
    ap = RingElement(round(a_vec),PQ);
    e_vec = randn([1,setup.N])*vpa(setup.sigma); % probably too small in practice
    ep = RingElement(round(e_vec),PQ);
    secret.q = PQ;
    bp = -ap*secret+ep+setup.P*secret^2;
    
    evk.c0 = bp;
    evk.c1 = ap;
    setup.evk = evk;
end