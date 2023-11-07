classdef (InferiorClasses = {?sym}) LWE_Ciphertext

    properties
        A {mustBeCoefficients(A)} = [];            % ciphertext tuple
        b (:,1) {mustBeCoefficients(b)} = [];      % ciphertext tuple
        N {mustBePositive,mustBeInteger}           % key dimension
        q {mustBeModulus(q)} = [];                 % modulus
    end

    methods
        % constructor
        function ct = LWE_Ciphertext(setup,A,b)
            arguments
                setup {mustBeA(setup,'GSW')}
                A
                b
            end
            ct.A = A;
            ct.b = b;
            ct.N = setup.N;
            ct.q = setup.q;
        end

        function ctnew = copyProperties(ct)
            ctnew = ct;
            ctnew.A = vpa([]);
            ctnew.b = vpa([]);
        end

        % decryption
        function me = Dec(ct,sk)
            arguments
                ct {mustBeA(ct,'LWE_Ciphertext')}
                sk (:,1) {mustBeCoefficients(sk)}
            end
            me = Mod(ct.b+ct.A*sk,ct.q);
        end

        % computation
        function ct = uminus(ct)
            ct.A = -ct.A;
            ct.b = -ct.b;
        end

        function result = plus(a,b)
            result = handleinputs(a,b,'plus');
            if isa(a,'LWE_Ciphertext') && isa(b,'LWE_Ciphertext')
                result.A = Mod(a.A+b.A,result.q);
                result.b = Mod(a.b+b.b,result.q);
            elseif isa(a,'LWE_Ciphertext')
                result.A = a.A;
                result.b = Mod(a.b+b,result.q);
            elseif isa(b,'LWE_Ciphertext')
                result.A = b.A;
                result.b = Mod(b.b+a,result.q);
            end
        end

        function result = minus(a,b)
            result = handleinputs(a,b,'plus');
            if isa(a,'LWE_Ciphertext') && isa(b,'LWE_Ciphertext')
                result.A = Mod(a.A-b.A,result.q);
                result.b = Mod(a.b-b.b,result.q);
            elseif isa(a,'LWE_Ciphertext')
                result.A = a.A;
                result.b = Mod(a.b-b,result.q);
            elseif isa(b,'LWE_Ciphertext')
                result.A = -b.A;
                result.b = Mod(a-b.b,result.q);
            end
        end

        function [result] = mtimes(a,b)
            result = handleinputs(a,b,'mult');            
            if isa(a,'LWE_Ciphertext')
                result.A = a.A*b;
                result.b = a.b*b;
            elseif isa(b,'LWE_Ciphertext')
                result.A = b.A*a;
                result.b = b.b*a;
            end     
        end
        
    end

    methods (Hidden)
        function c = handleinputs(a,b,operation)  

            % check inputs and prepare am empty ciphertext c
            a_isnumber = (isa(a,'sym') || isa(a,'numeric'));
            b_isnumber = (isa(b,'sym') || isa(b,'numeric'));
            if isa(a,'LWE_Ciphertext') && isa(b,'LWE_Ciphertext')
                if operation == "mult"
                    error('use GSW ciphertexts for multiplciation')
                end
                c = copyProperties(a);
                if a.N ~= b.N
                    error('dimension error: same dimensions expected')
                end
            elseif isa(a,'LWE_Ciphertext') && b_isnumber
                c = copyProperties(a);
                if sum(double(b-floor(b)),'all') ~= 0
                    error('invalid input: integer required')
                end
                if length(b) ~= 1
                    error('invalid input: scalar input required')
                end
            elseif a_isnumber && isa(b,'LWE_Ciphertext')
                c = copyProperties(b);
                if sum(double(a-floor(a)),'all') ~= 0
                    error('invalid input: integer required')
                end
                if length(a) ~= 1
                    error('invalid input: scalar input required')
                end
            end
        end
    end
end