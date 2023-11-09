classdef (InferiorClasses = {?sym}) Ciphertext

    properties
        c = [];                             % ciphertext
        q  {mustBeModulus(q)} = [];         % = modulus pk^2
        pk {mustBeModulus(pk)} = [];
    end

    methods
        % constructor
        function ct = Ciphertext(c,pk,q)
            arguments
                c
                pk
                q
            end
            ct.c  = c;
            ct.pk = pk;
            ct.q  = q;
        end

        function ctnew = copyProperties(ct)
            ctnew = ct;
            ctnew.c = vpa([]);
        end

        % computation
        function result = plus(a,b)
            result = handleinputs(a,b,'plus');
            if isa(a,'Ciphertext') && isa(b,'Ciphertext')
                result.c = mod(a.c.*b.c,result.q);
            elseif isa(a,'Ciphertext')
                tmp = Powermod(a.pk+1,b,result.q);
                result.c = mod(a.c.*tmp,result.q);
            elseif isa(b,'Ciphertext')
                tmp = Powermod(b.pk+1,a,result.q);
                result.c = mod(b.c.*tmp,result.q);
            end
        end

        function [result] = mtimes(a,b)
            [result,is_scalar] = handleinputs(a,b,'mult');            
            if is_scalar 
                if isa(a,'Ciphertext')
                    result.c = Powermod(a.c,b,a.q);
                elseif isa(b,'Ciphertext')
                    result.c = Powermod(b.c,a,b.q);
                end
            else % matrix multiplication
                if isa(a,'Ciphertext') 
                    result.c = vpa(ones(size(a.c,1),size(b,2)));
                    for i = 1:size(a.c,1)
                        for j = 1:size(b,2)
                            for k = 1:size(a.c,2)
                                tmp = Powermod(a.c(i,k),b(k,j),result.q);
                                result.c(i,j) = mod(result.c(i,j)*tmp,result.q);
                            end
                        end
                    end
                elseif isa(b,'Ciphertext')
                    result.c = vpa(ones(size(a,1),size(b.c,2)));
                    for i = 1:size(a,1)
                        for j = 1:size(b.c,2)
                            for k = 1:size(a,2)
                                tmp = Powermod(b.c(k,j),a(i,k),result.q);
                                result.c(i,j) = mod(result.c(i,j)*tmp,result.q);
                            end
                        end
                    end
                end
            end      
        end

        % concatenation
        function result = horzcat(a,b)
            result = copyProperties(a);
            result.c = [a.c, b.c];
        end

        function result = vertcat(a,b)
            result = copyProperties(a);
            result.c = [a.c;b.c];
        end
    end

    methods (Hidden)
        function [c,is_scalar] = handleinputs(a,b,operation)  

            % check inputs and prepare am empty ciphertext c
            a_isnumber = (isa(a,'sym') || isa(a,'numeric'));
            b_isnumber = (isa(b,'sym') || isa(b,'numeric'));
            if isa(a,'Ciphertext') && isa(b,'Ciphertext')
                if operation == "mult"
                    error('encrypted multiplication not supported')
                end
                c = copyProperties(a);
                a = a.c;
                b = b.c;
            elseif isa(a,'Ciphertext') && b_isnumber
                c = copyProperties(a);
                a = a.c;
                if sum(double(b-floor(b)),'all') ~= 0
                    error('invalid input: integer required')
                end
            elseif a_isnumber && isa(b,'Ciphertext')
                c = copyProperties(b);
                b = b.c;
                if sum(double(a-floor(a)),'all') ~= 0
                    error('invalid input: integer required')
                end
            end

            if operation == "plus"
                is_scalar = []; % not needed due to vectorization
                if any(size(a) ~= size(b))
                    error('dimension error: same dimensions expected')
                end
            end

            if operation == "mult"
                is_scalar = true;
                if any(size(a,2) ~= size(b,1))
                    error('dimension error: incompatible dimemnsions')
                end
                if size(a,2) ~= 1
                    is_scalar = false;
                end       
            end

        end
    end
end