classdef GSW
    % Bundles the cryptosystem configuration, provides key generation and
    % encryption utilities.
    properties
        N {mustBePositive,mustBeInteger}           % key dimension
        q {mustBeModulus(q)} = [];                 % modulus
        r {mustBePositive}                         % error bound
    end

    methods
        % constructor
        function setup = GSW(N,q,r)
            arguments
                N
                q
                r
            end
            setup.N = N;
            setup.q = q;
            setup.r = r;
        end

        function sk = KeyGen(setup)
            % Generates the secret key
            arguments
                setup {mustBeA(setup,'GSW')}
            end
            % NOTE: The randomness and emulated distributions do not satisfy
            %       cryptographic requirements (see comments below)!
            sk = Mod(round(rand([setup.N, 1])*setup.q),setup.q);
        end

        function ct = Enc(setup,m,sk,options)
            arguments
                setup    {mustBeA(setup,'GSW')}
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

        function [ct,ct0,ct_tuple] = Enc2(setup,m,sk,options)
            arguments
                setup    {mustBeA(setup,'GSW')}
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
                        ct0 = Enc(setup,m0,sk,'error',false);
                    else
                        ct0 = Enc(setup,m0,sk);
                    end
                    ct_tuple = Mod(m*R+[ct0.b ct0.A],setup.q);
                end
            else 
                %ct_tuple = vpa(zeros(double(log10(setup.q))*...
                %           (setup.N+1), setup.N+1, n1, n2)); % initialize
                for i = 1:n1
                    for j = 1:n2
                        m0 = vpa(zeros(double(lq*(setup.N+1)),1));
                        ct0 = Enc(setup,m0,sk);
                        ct_tuple{i,j} = Mod(m(i,j)*R ...
                            +[ct0.b ct0.A],setup.q);
                    end
                end
            end
            ct = GSW_Ciphertext(setup,ct_tuple);
        end

    end

end