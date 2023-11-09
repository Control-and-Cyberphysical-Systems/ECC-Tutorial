function [ct,ct0,ct_tuple] = GSW(setup,m,sk,options)
% GSW encryption
    arguments
        setup    {mustBeA(setup,'Setup')}
        m        {mustBeCoefficients(m)}
        sk (:,1) {mustBeCoefficients(sk)}
        options.debug = false;
        options.error = true;
    end
    
    [n1,n2] = size(m);
    lq = log10(setup.q);
    base = vpa(10);
    R = kron(power(base,(0:1:lq-1)'),eye(setup.N+1));

    if n1 == 1 && n2 == 1 % scalar case
        if options.debug == true
            ct_tuple = m*R;
        else
            m0 = vpa(zeros(double(lq*(setup.N+1)),1));
            if options.error == false
                ct0 = LWE(setup,m0,sk,'error',false);
            else
                ct0 = LWE(setup,m0,sk);
            end
            ct_tuple = Mod(m*R+[ct0.b ct0.A],setup.q);
        end
    else 
        %ct_tuple = vpa(zeros(double(log10(setup.q))*...
        %           (setup.N+1), setup.N+1, n1, n2)); % initialize
        for i = 1:n1
            for j = 1:n2
                m0 = vpa(zeros(double(lq*(setup.N+1)),1));
                ct0 = LWE(setup,m0,sk);
                ct_tuple{i,j} = Mod(m(i,j)*R ...
                    +[ct0.b ct0.A],setup.q);
            end
        end
    end
    ct = GSW_Ciphertext(setup,ct_tuple);
end