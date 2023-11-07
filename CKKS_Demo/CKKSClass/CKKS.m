classdef CKKS
    % Bundles the cryptosystem configuration, provides key generation and
    % encryption utilities.
    properties
        N {mustBePositive,mustBeInteger}           % ring dimension
        q0 {mustBeModulus(q0)} = [];               % base modulus
        s  {mustBeScalingFactor(s)} = [];          % scaling factor
        L {mustBePositive,mustBeInteger} = [];     % ciphertext levels
        P {mustBePositive} = [];                   % huge scaling factor
        sigma {mustBeGreaterThanOrEqual(sigma,0)}  % standard deviation

        evk                                        % stored here for convenient * overloading
    end

    properties (Dependent)
        Q                                          % starting modulus
    end

    methods
        % constructor
        function setup = CKKS(N,q0,s,L,P,options)
            % Inputs:
            %   N: integer
            %       number of coefficients in polynomials
            %   q0: integer
            %       base modulus
            %   s:
            %       scaling factor
            %   L: integer
            %       maximum number of levels
            %   P:
            %       large number related to the evaluation key
            %   options: struct
            %       sigma:
            %           standard deviation for Gaussian sampling (default: 3.19)
            arguments
                N
                q0
                s
                L
                P
                options.sigma = 3.19
            end
            setup.N = N;
            setup.q0 = q0;
            setup.s = s;
            setup.L = L;
            setup.P = P;
            setup.sigma = options.sigma;
        end

        function Q = get.Q(setup)
            Q = ceil(setup.q0*setup.s^setup.L);
        end

        function [pk,sk,setup] = KeyGen(setup)
            % Generate the public, secret and evaluation key. The latter is
            % stored as a property of the returned `CKKS` instance `setup`.
            arguments
                setup {mustBeA(setup,'CKKS')}
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

        function swk = swichtingKey(setup,sk,skprime)
            % Generate the key that enables switching ciphertexts that are valid
            % for decryption under `sk` to become valid for decryption under
            % `skprime` (decoupled from KeyGen).
            arguments
                setup {mustBeA(setup,'CKKS')}
                sk {mustBeA(sk,'RingElement')}
                skprime {mustBeA(skprime,'RingElement')}
            end
            PQ = setup.P*setup.Q;
            skprime.q = PQ;
            sk.q = PQ;

            a_vec = vpa(rand([1,setup.N]))*PQ-PQ/2;
            as = RingElement(round(a_vec),PQ);
            e_vec = randn([1,setup.N])*vpa(setup.sigma); % probably too small in practice
            es = RingElement(round(e_vec),PQ);
            bs = -as*skprime+es+setup.P*sk;

            swk.c0 = bs;
            swk.c1 = as;
        end

        function ct = Enc(setup,m,pk)
            arguments
                setup {mustBeA(setup,'CKKS')}
                m {mustBeA(m,'Encoding')}
                pk {mustBeKey(pk)}
            end
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
            ct = Ciphertext(setup,m,c0,c1);
        end

    end

end