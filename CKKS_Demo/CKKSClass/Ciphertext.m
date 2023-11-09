classdef (InferiorClasses = {?sym,?RingElement}) Ciphertext
    % Represents a CKKS ciphertext
    properties
        c0 {mustBeCtTuple(c0)}                     % ciphertext tuple
        c1 {mustBeCtTuple(c1)}

        level {mustBeInteger} = [];                % current level
        s  {mustBeScalingFactor(s)} = [];          % scaling factor
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

    properties (Hidden)
        type                                    % used for decoding
        is_scalar                               % used for distinguishing
    end

    methods
        % constructor
        function ct = Ciphertext(setup,m,c0,c1)
            % Construct a new ciphertext `(c0, c1)` that encrypts the message
            % `m`
            arguments
                setup {mustBeA(setup,'Setup')}
                m
                c0 = [];
                c1 = [];
            end
            % setup properties
            ct.q0 = setup.q0;
            ct.P = setup.P;
            ct.L = setup.L;
            ct.Q = setup.Q;
            ct.evk = setup.evk; % # sneaky
            ct.s = setup.s;
            ct.level = setup.L;
            
            % encoding
            ct.type = m.type; 
            if ~iscell(c0) % is scalar?
                ct.is_scalar = 1;
            else
                ct.is_scalar = 0;
            end

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

        function ct = set.c0(ct,c0)
            if ~iscell(c0) % is scalar?
                ct.is_scalar = 1;
            else
                ct.is_scalar = 0;
            end
            ct.c0 = c0;
        end

        function ct = set.q(ct,q)
            if ct.is_scalar
                ct.c0.q = q;
            else
                ct.c0.q{1,1} = q;
            end
        end

        function q = get.q(ct)
            if ct.is_scalar
                q = ct.c0.q;
            else
                q = ct.c0{1,1}.q;
            end
        end

        function N = get.N(ct)
            % The ring dimension of the contained ring elements
            if ct.is_scalar
                N = ct.c0.N;
            else
                N = ct.c0{1,1}.N;
            end
        end   
        
        % computations
        function ct = uminus(ct)
            if ct.is_scalar
                ct.c0 = -ct.c0;
                ct.c1 = -ct.c1;
            else
                for i = 1:size(ct.c0,1)
                    for j = 1:size(ct.c0,2)
                        ct.c0{i,j} = -ct.c0{i,j};
                        ct.c1{i,j} = -ct.c1{i,j};
                    end
                end
            end
        end

        function c = plus(a,b)
            % Encrypted addition between two ciphertexts or a ciphertext 
            % and an unencrypted RingElement
            [a,b,c] = handleinputs(a,b,'constant');
            if isa(a,'Ciphertext') && isa(b,'Ciphertext')
                if a.is_scalar == 1 && b.is_scalar == 1
                    c.c0 = a.c0+b.c0;
                    c.c1 = a.c1+b.c1;
                elseif ~a.is_scalar && ~b.is_scalar
                    if any(size(a.c0)-size(b.c0))
                       error('incompatible dimensions')
                    end
                    c.c0 = {};
                    c.c1 = {};
                    c.c0{1,1} = RingElement(vpa(0),a.q);
                    c.c1{1,1} = RingElement(vpa(0),a.q);
                    for i = 1:size(a.c0,1)
                        for j = 1:size(a.c0,2)
                            c.c0{i,j} = a.c0{i,j}+b.c0{i,j};
                            c.c1{i,j} = a.c1{i,j}+b.c1{i,j};
                        end
                    end
                else
                    error('incompatible dimensions')
                end
            elseif isa(a,'Ciphertext') && isa(b,'RingElement')
                c.c0 = a.c0+b;
                c.c1 = a.c1;
            elseif isa(a,'RingElement') && isa(b,'Ciphertext')
                c.c0 = b.c0+a;
                c.c1 = b.c1;
            else
                error('this case is not implemented')
            end
        end

        function c = minus(a,b)
            [a,b,c] = handleinputs(a,b,'constant');
            if isa(a,'Ciphertext') && isa(b,'Ciphertext')
                if a.is_scalar == 1 && b.is_scalar == 1
                    c.c0 = a.c0-b.c0;
                    c.c1 = a.c1-b.c1;
                elseif ~a.is_scalar && ~b.is_scalar
                    if any(size(a.c0)-size(b.c0))
                       error('incompatible dimensions')
                    end
                    c.c0 = {};
                    c.c1 = {};
                    c.c0{1,1} = RingElement(vpa(0),a.q);
                    c.c1{1,1} = RingElement(vpa(0),a.q);
                    for i = 1:size(a.c0,1)
                        for j = 1:size(a.c0,2)
                            c.c0{i,j} = a.c0{i,j}-b.c0{i,j};
                            c.c1{i,j} = a.c1{i,j}-b.c1{i,j};
                        end
                    end
                else
                    error('incompatible dimensions')
                end
            elseif isa(a,'Ciphertext') && isa(b,'RingElement')
                c.c0 = a.c0-b;
                c.c1 = a.c1;
            elseif isa(a,'RingElement') && isa(b,'Ciphertext')
                c.c0 = b.c0-a;
                c.c1 = b.c1;
            else
                error('this case is not implemented')
            end
        end

        function c = mtimes(a,b)
            % Encrypted muliplication between two ciphertexts or a ciphertext
            % and an unencrypted RingElement
            [a,b,c] = handleinputs(a,b,'increase');
            if isa(a,'Ciphertext') && isa(b,'Ciphertext')
                if a.is_scalar == 1 && b.is_scalar == 1
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

                else
                    c = matmult(a,b,c);
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
            else
                error('this case is not implemented')
            end
        end

        function c = matmult(a,b,c) % naive implementation to avoid rotaions
            if a.type == "packed"
                error('expected single scalar based encoding')
            end
            if size(a.c0,2)~=size(b.c0,1)
                error('incompatible dimesions')
            end
            c.c0={};
            c.c1={};
            for i = 1:size(a.c0,1)
                for k = 1:size(b.c0,2)
                    flag = 1;
                    for j = 1:size(a.c0,2)
                        a_ctij = a;                 % copy
                        a_ctij.c0 = a_ctij.c0{i,j}; % delete all but ij ct
                        a_ctij.c1 = a_ctij.c1{i,j};
                        b_ctjk = b;
                        b_ctjk.c0 = b_ctjk.c0{j,k};  
                        b_ctjk.c1 = b_ctjk.c1{j,k};

                        ctmp = a_ctij*b_ctjk;

                        if flag == 1
                            c.c0{i,k} = RingElement(vpa(0),ctmp.q);
                            c.c1{i,k} = RingElement(vpa(0),ctmp.q);
                            flag = 0;
                            c.s = ctmp.s;
                            c.level = ctmp.level;
                        end

                        c.c0{i,k} = c.c0{i,k}+ctmp.c0;
                        c.c1{i,k} = c.c1{i,k}+ctmp.c1;
                    end
                end
            end
            if size(c.c0) == 1
                c.c0 = c.c0{1};
                c.c1 = c.c1{1};
            end
        end

        function [c,powct] = mpower(ct,b)
            % (faster) Encrypted exponentiation of a ciphertext `ct` by the power `b`
            arguments
                ct {mustBeA(ct,'Ciphertext')}
                b {mustBeScalarOrEmpty}
            end
            if ~ct.is_scalar
                error('keyswitching not implemented for non-scalar ciphertexts');
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
            if ct.is_scalar
            ct.c0 = floor(ct.c0);
            ct.c1 = floor(ct.c1);
            else
                for i = 1:size(ct.c0,1)
                    for j = 1:size(ct.c0,2)
                        ct.c0{i,j} = floor(ct.c0{i,j});
                        ct.c1{i,j} = floor(ct.c1{i,j});
                    end
                end
            end
        end

        function ct = ceil(ct)
            if ct.is_scalar
            ct.c0 = ceil(ct.c0);
            ct.c1 = ceil(ct.c1);
            else
                for i = 1:size(ct.c0,1)
                    for j = 1:size(ct.c0,2)
                        ct.c0{i,j} = ceil(ct.c0{i,j});
                        ct.c1{i,j} = ceil(ct.c1{i,j});
                    end
                end
            end
        end

        function ct = round(ct)
            if ct.is_scalar
            ct.c0 = round(ct.c0);
            ct.c1 = round(ct.c1);
            else
                for i = 1:size(ct.c0,1)
                    for j = 1:size(ct.c0,2)
                        ct.c0{i,j} = round(ct.c0{i,j});
                        ct.c1{i,j} = round(ct.c1{i,j});
                    end
                end
            end
        end

        % modulus/scale
        function ct = modreduce(ct,qprime)
            arguments
                ct {mustBeA(ct,'Ciphertext')}
                qprime {mustBeModulus(qprime)}
            end
            if ct.is_scalar
                ct.c0 = modreduce(ct.c0,qprime);
                ct.c1 = modreduce(ct.c1,qprime);
            else
                for i = 1:size(ct.c0,1)
                    for j = 1:size(ct.c0,2)
                        ct.c0{i,j} = modreduce(ct.c0{i,j},qprime);
                        ct.c1{i,j} = modreduce(ct.c1{i,j},qprime);
                    end
                end
            end
            if double(ct.q<ct.q0)
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

            if ct.is_scalar
                ct.c0 = rescale(ct.c0,b);
                ct.c1 = rescale(ct.c1,b);
            else
                for i = 1:size(ct.c0,1)
                    for j = 1:size(ct.c0,2)
                        ct.c0{i,j} = rescale(ct.c0{i,j},b);
                        ct.c1{i,j} = rescale(ct.c1{i,j},b);
                    end
                end
            end
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
            if ~ct.is_scalar
                error('keyswitching not implemented for non-scalar ciphertexts');
            end
            ct.c1.q = ct.P*ct.Q;
            ks.c0 = round(1/ct.P*ct.c1*swk.c0);
            ks.c1 = round(1/ct.P*ct.c1*swk.c1);
            ks.c0 = modreduce(ks.c0,ct.q);
            ks.c1 = modreduce(ks.c1,ct.q);

            ct.c0 = ct.c0+ks.c0;
            ct.c1 = ks.c1;
        end

        function [a,b] = extractcoef(ct)
            % Extract the coefficients of the contained ring elements
            if ~ct.is_scalar
                error('keyswitching not implemented for non-scalar ciphertexts');
            end
            b = ct.c0.coefficients(1);
            a = fliplr(ct.c1.coefficients);
        end

        function ct_transposed = ctranspose(ct)
            if ~ct.is_scalar
                ct_transposed = ct;
                ct_transposed.c0 = {};
                ct_transposed.c1 = {};
                ct_transposed.c0{1,1} = RingElement(0,ct.q);
                ct_transposed.c1{1,1} = RingElement(0,ct.q);
                for i = 1:size(ct.c0,1)
                    for j = 1:size(ct.c1,2)
                        ct_transposed.c0{j,i} = ct.c0{i,j};
                        ct_transposed.c1{j,i} = ct.c1{i,j};
                    end
                end
            end
        end

        function ct = vertcat(ct,ct2)
            if ct.level ~= ct2.level || ct.q ~= ct2.q || ct.s ~= ct2.s ...
               || ct.N ~= ct2.N
                error("Cipherexts do not share properties." + ...
                    " Concatenation not possible")
            end
            ct.c0 = [ct.c0; ct2.c0];
            ct.c1 = [ct.c1; ct2.c1];
        end
        
        function ct = horzcat(ct,ct2)
            if ct.level ~= ct2.level || ct.q ~= ct2.q || ct.s ~= ct2.s ...
               || ct.N ~= ct2.N
                error("Cipherexts do not share properties." + ...
                    " Concatenation not possible")
            end
            ct.c0 = [ct.c0 ct2.c0];
            ct.c1 = [ct.c1 ct2.c1];
        end

        % NOTE: The following are not implemented here
        %   - rotation
        %   - automorphism + keyswitch

    end

    methods (Hidden)
        function [a,b,c] = handleinputs(a,b,scaling)
            if isa(a,'Ciphertext') && isa(b,'Ciphertext') 
               mode = a.adapt;
               c = copyProperties(a); % retain setup properties
               c.level = min(a.level,b.level);
            elseif isa(a,'Ciphertext') && isa(b,'RingElement')
               c = copyProperties(a); % retain setup properties
               mode = a.adapt;
               if isempty(b.s)
                   error('Multiplication with non-encoded elements undefined')
               end
            elseif isa(a,'RingElement') && isa(b,'Ciphertext')
               mode = b.adapt;
               c = copyProperties(b);
               if isempty(a.s)
                   error('Multiplication with non-encoded elements undefined')
               end
            else
                error("invalid input")
            end

            [a,b] = adaptModScale(a,b,mode,scaling);
            % compute the scaling
            if strcmp(scaling,'increase')
                c.s = a.s*b.s;
            elseif strcmp(scaling,'constant')
                c.s = a.s;
            end
            c.q = a.q; % set the the modulus

        end
    end

end