function ct = LWE(setup,m,sk,options)
% LWE encryption
    arguments
        setup    {mustBeA(setup,'Setup')}
        m  (:,1) {mustBeCoefficients(m)}
        sk (:,1) {mustBeCoefficients(sk)}
        options.debug = false;
        options.error = true;
    end
    
    n = length(m);
    
    if options.debug 
        A = zeros([n, setup.N]);
        e = zeros([n,1]);
    else
        A  = Mod(round(rand([n, setup.N])*setup.q),setup.q);
        if options.error == false
            e = vpa(zeros(n,1));
        else
            e  = round(rand([n,1])*setup.r-setup.r/2);
        end
    end
    b  = Mod(-A*sk + m + e,setup.q);

    ct = LWE_Ciphertext(setup,A,b);
end