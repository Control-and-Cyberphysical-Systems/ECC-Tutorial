classdef (InferiorClasses = {?sym,?RingElement,?Encoding}) Ciphertext
    % Represents a CKKS ciphertext. Notably, this consists of two ring elements.
    % Provides encrypted operations, key switching and decryption.
    properties
        c0 RingElement = RingElement();            % ciphertext tuple
        c1 RingElement = RingElement();

        level {mustBeInteger} = [];                % current level
        s  {mustBeScalingFactor(s)} = [];          % scaling factor
        type {mustBeEncodingType(type)} = [];      % encoding type
        adapt {mustBeAdaptionMode(adapt)} = 'auto' % modulo and scaling mode
    end

    properties (SetAccess = immutable)
        q0 {mustBeModulus(q0)} = [];            % base modulus
        P {mustBePositive} = [];                % huge scaling
        L {mustBePositive,mustBeInteger} = [];  % total levels
        Q {mustBeModulus(Q)}                    % starting modulus
        evk                                     % for convenient * overloading
    end

    properties (Dependent)
        q  % current modulus
        N  % ring dimension
    end

    methods
        % constructor
        function ct = Ciphertext(setup,m,c0,c1)
            % Construct a new ciphertext `(c0, c1)` that encrypts the message
            % `m`. Note that `m` (a Encoding) is only used to extract the
            % encoding properties (scaling and encoding type).
            arguments
                setup {mustBeA(setup,'CKKS')}
                m {mustBeA(m,'Encoding')}
                c0 {mustBeA(c0,'RingElement')}
                c1 {mustBeA(c1,'RingElement')}
            end
            % setup properties
            ct.q0 = setup.q0;
            ct.P = setup.P;
            ct.L = setup.L;
            ct.Q = setup.Q;
            ct.evk = setup.evk; % # sneaky

            % encoding properties
            ct.level = ct.L;
            ct.s = m.s;
            ct.type = m.type;

            % ciphertext tuple
            ct.c0 = c0;
            ct.c1 = c1;
        end

        function ctnew = copyProperties(ct)
            % Generate an empty ciphertext with the otherwise same properties as
            % this one.
            ctnew = ct;
            ctnew.c0 = RingElement([],ct.q);
            ctnew.c1 = RingElement([],ct.q);
        end

        function q = get.q(ct)
            % The modulus of the contained ring elements
            q = ct.c0.q;
        end

        function N = get.N(ct)
            % The ring dimension of the contained ring elements
            N = ct.c0.N;
        end

        % decryption
        function me = Dec(ct,sk)
            % Decrypt using the secret key `sk`. Returns a Encoding.
            arguments
                ct {mustBeA(ct,'Ciphertext')}
                sk {mustBeA(sk,'RingElement')}
            end
            sk.q = ct.q;    % ternary (no reduction needed)
            me_ring = ct.c0+ct.c1*sk;
            me = Encoding(me_ring.coefficients,ct.q,ct.s,ct.type,ct.adapt);
        end

        % computation
        function ct = uminus(ct)
            ct.c0 = -ct.c0;
            ct.c1 = -ct.c1;
        end

        function c = plus(a,b)
            % Encrypted addition between two ciphertexts or a ciphertext and an
            % unencrypted RingElement
            [a,b,c] = handleinputs(a,b,'constant');
            if isa(a,'Ciphertext') && isa(b,'Ciphertext')
                c.c0 = a.c0+b.c0;
                c.c1 = a.c1+b.c1;
            elseif isa(a,'Ciphertext') && isa(b,'RingElement')
                c.c0 = a.c0+b;
                c.c1 = a.c1;
            elseif isa(a,'RingElement') && isa(b,'Ciphertext')
                c.c0 = b.c0+a;
                c.c1 = b.c1;
            end
        end

        function c = minus(a,b)
            % Encrypted subtraction between two ciphertexts or a ciphertext and
            % an unencrypted RingElement
            [a,b,c] = handleinputs(a,b,'constant');
            if isa(a,'Ciphertext') && isa(b,'Ciphertext')
                c.c0 = a.c0-b.c0;
                c.c1 = a.c1-b.c1;
            elseif isa(a,'Ciphertext') && isa(b,'RingElement')
                c.c0 = a.c0-b;
                c.c1 = a.c1;
            elseif isa(a,'RingElement') && isa(b,'Ciphertext')
                c.c0 = b.c0-a;
                c.c1 = b.c1;
            end
        end

        function c = mtimes(a,b)
            % Encrypted muliplication between two ciphertexts or a ciphertext
            % and an unencrypted RingElement
            [a,b,c] = handleinputs(a,b,'increase');
            if isa(a,'Ciphertext') && isa(b,'Ciphertext')
                % multiplication
                d0 = a.c0*b.c0;
                d1 = a.c0*b.c1+a.c1*b.c0;
                d2 = a.c1*b.c1;

                % relinearization
                d2.q = c.P*c.Q;
                d2prime.c0 = round(1/c.P*d2*c.evk.c0);
                d2prime.c1 = round(1/c.P*d2*c.evk.c1);
                d2prime.c0 = modreduce(d2prime.c0,c.q);
                d2prime.c1 = modreduce(d2prime.c1,c.q);

                c.c0 = d0+d2prime.c0;
                c.c1 = d1+d2prime.c1;

                if strcmp(c.adapt,'auto')
                    c = rescale(c,a.s);
                end
            elseif isa(a,'Ciphertext') && isa(b,'RingElement')
                c.c0 = a.c0*b;
                c.c1 = a.c1*b;
                if strcmp(c.adapt,'auto')
                    c = rescale(c,a.s);
                end
            elseif isa(a,'RingElement') && isa(b,'Ciphertext')
                c.c0 = b.c0*a;
                c.c1 = b.c1*a;
                if strcmp(c.adapt,'auto')
                    c = rescale(c,b.s);
                end
            end
        end

        function [c,powct] = mpower(ct,b)
            % (faster) Encrypted exponentiation of a ciphertext `ct` by the power `b`
            arguments
                ct {mustBeA(ct,'Ciphertext')}
                b {mustBeScalarOrEmpty}
            end
            bin_exp = floor(log2(b));
            powct = cell(bin_exp,1);
            powct{1} = ct;
            for i = 1:bin_exp
                powct{i+1} = powct{i}*powct{i};
            end
            selection = de2bi(b-2^bin_exp,bin_exp);
            c = powct{end};
            for i = 1:length(selection)
                if selection(i) == 1
                    c = c*powct{i};
                end
            end
        end

        % rounding etc.
        function ct = floor(ct)
            ct.c0 = ceil(ct.c0);
            ct.c1 = ceil(ct.c1);
        end

        function ct = ceil(ct)
            ct.c0 = ceil(ct.c0);
            ct.c1 = ceil(ct.c1);
        end

        function ct = round(ct)
            ct.c0 = round(ct.c0);
            ct.c1 = round(ct.c1);
        end

        % modulus/scale
        function ct = modreduce(ct,qprime)
            % modreduce for each contained ring element
            arguments
                ct {mustBeA(ct,'Ciphertext')}
                qprime {mustBeModulus(qprime)}
            end
            ct.c0 = modreduce(ct.c0,qprime);
            ct.c1 = modreduce(ct.c1,qprime);
            if ct.q <= ct.q0
                warning("current modulus is below the base modulus." + ...
                    " Decryption may be incorrect")
            end
        end

        function ct = rescale(ct,b)
            % Rescale each contained ring element. This also reduces the level
            % by one.
            arguments
                ct {mustBeA(ct,'Ciphertext')}
                b {mustBePositive,mustBeScalarOrEmpty}
            end
            ct.c0 = rescale(ct.c0,b);
            ct.c1 = rescale(ct.c1,b);
            ct.s = round(ct.s/b);
            if double(b ~= 1)
                ct.level = ct.level - 1;
            end
            if ct.level < 0
                warning("ciphertext level is below 0. Decryption may" +...
                    " be incorrect")
            end
        end

        % keys
        function ct = keyswitch(ct,swk)
            % Return a new ciphertext encrypting (conceptually) the same message
            % which is valid for decryption under a different secret key (see
            % `CKKS.switchingKey` for details).
            arguments
                ct {mustBeA(ct,'Ciphertext')}
                swk {mustBeKey(swk)}
            end
            ct.c1.q = ct.P*ct.Q;
            ks.c0 = round(1/ct.P*ct.c1*swk.c0);
            ks.c1 = round(1/ct.P*ct.c1*swk.c1);
            ks.c0 = modreduce(ks.c0,ct.c0.q);
            ks.c1 = modreduce(ks.c1,ct.c0.q);

            ct.c0 = ct.c0+ks.c0;
            ct.c1 = ks.c1;
        end

        function [a,b] = extractcoef(ct)
            % Extract the coefficients of the contained ring elements
            b = ct.c0.coefficients(1);
            a = fliplr(ct.c1.coefficients);
        end

        % NOTE: The following are not implemented here
        %   - rotation
        %   - automorphism + keyswitch

    end

    methods (Hidden)
        function [a,b,c] = handleinputs(a,b,scaling)
            if isa(a,'Ciphertext') && isa(b,'Ciphertext') 
               c = copyProperties(a); % retain setup properties
               c.level = min(a.level,b.level);
            elseif isa(a,'Ciphertext') && isa(b,'Encoding')
               c = copyProperties(a); % retain setup properties
            elseif isa(a,'Encoding') && isa(b,'Ciphertext')
               c = copyProperties(b);
            else
                error('invalid input: ciphertexts and encodings expected')
            end
            [a,b,c.adapt] = adaptModScale(a,b,scaling);
            % compute the scaling
            if strcmp(scaling,'increase')
                c.s = a.s*b.s;
            elseif strcmp(scaling,'constant')
                c.s = a.s;
            end
            c.c0.q = a.q;
            % initialize for computation
            if isa(a,'Ciphertext') && isa(b,'Encoding')
                b = RingElement(b.coefficients,b.q);
            elseif isa(a,'Encoding') && isa(b,'Ciphertext')
                a = RingElement(a.coefficients,a.q);
            end
        end
    end

end