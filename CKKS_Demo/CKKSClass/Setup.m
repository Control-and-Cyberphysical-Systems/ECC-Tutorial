classdef Setup
    % Bundles the cryptosystem configuration
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
        function setup = Setup(N,q0,s,L,P,options)
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
                P = q0*q0;
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

        function swk = swichtingKey(setup,sk,skprime)
            % Generate the key that enables switching ciphertexts that are valid
            % for decryption under `sk` to become valid for decryption under
            % `skprime` (decoupled from KeyGen).
            arguments
                setup {mustBeA(setup,'Setup')}
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

    end
end