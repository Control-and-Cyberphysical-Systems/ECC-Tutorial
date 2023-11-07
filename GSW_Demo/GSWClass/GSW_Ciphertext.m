classdef (InferiorClasses = {?sym,?LWE_Ciphertext}) GSW_Ciphertext

    properties
        c {mustBeCoefficients(c)} = [];            % ciphertext
        N {mustBePositive,mustBeInteger}           % key dimension
        q {mustBeModulus(q)} = [];                 % modulus
    end

    methods
        % constructor
        function ct = GSW_Ciphertext(setup,c)
            arguments
                setup {mustBeA(setup,'GSW')}
                c
            end
            ct.c = c;
            ct.N = setup.N;
            ct.q = setup.q;
        end

        function ctnew = copyProperties(ct)
            ctnew = ct;
            ctnew.c = vpa([]);
        end

        % decryption
        function me = Dec(ct,sk)
            arguments
                ct {mustBeA(ct,'GSW_Ciphertext')}
                sk
            end
            % select one element  and decrypt
            me = Mod(ct.c(1,:)*[1;sk],ct.q);
        end

        % computation
        function ct = uminus(ct)
            ct.c = -ct.c;
        end

        function result = plus(a,b)
            result = handleinputs(a,b,'plus');

            if isa(a,'GSW_Ciphertext') && isa(b,'GSW_Ciphertext')
                result.c = Mod(a.c+b.c,result.q);
            elseif isa(a,'GSW_Ciphertext')
                lq = log10(a.q);
                base = vpa(10);
                R = kron(power(base,(0:1:lq-1)'),eye(a.N+1));
                result.c = Mod(a.c+R*b,result.q);
            elseif isa(b,'GSW_Ciphertext')
                lq = log10(b.q);
                base = vpa(10);
                R = kron(power(base,(0:1:lq-1)'),eye(b.N+1));
                result.c = Mod(b.c+R*a,result.q);
            end
        end

        function result = minus(a,b)
            result = handleinputs(a,b,'plus');
            
            if isa(a,'GSW_Ciphertext') && isa(b,'GSW_Ciphertext')
                result.c = Mod(a.c-b.c,result.q);
            elseif isa(a,'GSW_Ciphertext')
                lq = log10(a.q);
                base = vpa(10);
                R = kron(power(base,(0:1:lq-1)'),eye(a.N+1));
                result.c = Mod(a.c-R*b,result.q);
            elseif isa(b,'GSW_Ciphertext')
                lq = log10(b.q);
                base = vpa(10);
                R = kron(power(base,(0:1:lq-1)'),eye(b.N+1));
                result.c = Mod(-b.c+R*a,result.q);
            end
        end

        function result = mtimes(a,b)
            [result,is_scalar] = handleinputs(a,b,'mult');
            
            if is_scalar
                if isa(a,'GSW_Ciphertext') && isa(b,'LWE_Ciphertext')
                    tmp = Mod(Decomp(b)*a.c,result.q);
                    result.b = tmp(:,1);
                    result.A = tmp(:,2:end);
                elseif isa(a,'LWE_Ciphertext') && isa(b,'GSW_Ciphertext')
                    tmp = Mod(Decomp(a)*b.c,result.q);
                    result.b = tmp(:,1);
                    result.A = tmp(:,2:end);
                elseif isa(a,'GSW_Ciphertext')
                    result.c = Mod(a.c*b,result.q);
                elseif isa(b,'GSW_Ciphertext')
                    result.c = Mod(b.c*a,result.q);
                end
            end

            if ~is_scalar % matrix multiplication
                if isa(a,'GSW_Ciphertext') && isa(b,'LWE_Ciphertext')
                    Db = Decomp(b);
                    if size(Db,2) ~= size(a.c{1,1},1) || ....
                       size(Db,1) ~= size(a.c,1)
                        error('dimension error: ciphertexts not compatible')
                    end
                    for i = 1:size(a.c,1)
                        result.b(i) = vpa(0);
                        result.A(i,:) = vpa(zeros(1,b.N));
                        for j = 1:size(a.c,2)
                            tmp = Mod(Db(j,:)*a.c{i,j},result.q);
                            result.b(i) = Mod(result.b(i)+tmp(1),result.q);
                            result.A(i,:) = Mod(result.A(i,:)+tmp(2:end),result.q);
                        end

                    end
                else
                    error('unexpected input')
                end
            end

        end

        function ct_LWE = extractLWE(ct_GSW)
            tmp = ct_GSW.c(1,:);
            b = tmp(:,1);
            A = tmp(:,2:end);
            setup = GSW(ct_GSW.N,ct_GSW.q,vpa(1));
            ct_LWE = LWE_Ciphertext(setup,A,b);
        end

    end

    methods (Hidden)
        function [c,is_scalar] = handleinputs(a,b,operation)

            % check inputs and prepare am empty ciphertext c
            a_isnumber = (isa(a,'sym') || isa(a,'numeric'));
            b_isnumber = (isa(b,'sym') || isa(b,'numeric'));
            is_scalar = true;

            if operation == "plus"
                if isa(a,'GSW_Ciphertext') && isa(b,'GSW_Ciphertext')
                    if any(size(a.c) ~= size(b.c))
                        error('dimension error: same dimensions expected')
                    end
                    c = copyProperties(a);
                elseif isa(a,'GSW_Ciphertext') && isa(b,'LWE_Ciphertext') ...
                        || isa(b,'LWE_Ciphertext') && isa(b,'GSW_Ciphertext')
                    error('Addition between GSW and LWE ciphertexts not supported')
                elseif isa(a,'GSW_Ciphertext') && b_isnumber
                    c = copyProperties(a);
                    if sum(double(b-floor(b)),'all') ~= 0
                        error('invalid input: integer required')
                    end
                    if length(b) ~= 1
                        error('dimension error: scalar expected')
                    end
                elseif a_isnumber && isa(b,'GSW_Ciphertext')
                    c = copyProperties(b);
                    if sum(double(a-floor(a)),'all') ~= 0
                        error('invalid input: integer required')
                    end
                    if length(a) ~=1
                        error('dimension error: scalar expected')
                    end
                else
                    error(['Type error: GSW ciphertext can only be' ...
                        'added to GSW ciphertexts or scalars'])
                end
            end

            if operation == "mult"
                if isa(a,'GSW_Ciphertext') && isa(b,'GSW_Ciphertext')
                    error('ciphertext multiplication type not supported \n%s', ...
                        'GSW * GSW -> GSW (not implemented) \n%s', ...
                        'Decomp(LWE) * GSW -> LWE (here)')
                elseif isa(a,'LWE_Ciphertext') && isa(b,'GSW_Ciphertext')
                    c = copyProperties(a);
                    is_scalar = ~isa(b.c,'cell');
                elseif isa(a,'GSW_Ciphertext') && isa(b,'LWE_Ciphertext')
                    c = copyProperties(b);
                    is_scalar = ~isa(a.c,'cell');
                elseif isa(a,'GSW_Ciphertext') && b_isnumber
                    c = copyProperties(a);
                    is_scalar = ~isa(a.c,'cell');
                    if sum(double(b-floor(b)),'all') ~= 0
                        error('invalid input: integer required')
                    end
                    if length(b) ~= 1
                        error('dimension error: scalar expected')
                    end
                elseif a_isnumber && isa(b,'GSW_Ciphertext')
                    c = copyProperties(b);
                    is_scalar = ~isa(b.c,'cell');
                    if sum(double(a-floor(a)),'all') ~= 0
                        error('invalid input: integer required')
                    end
                    if length(a) ~=1
                        error('dimension error: scalar expected')
                    end
                end
            end

        end
    end
end