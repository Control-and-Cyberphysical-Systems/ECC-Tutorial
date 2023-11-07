classdef (InferiorClasses = {?sym}) RingElement
    % Representation of elements from polynomial rings with optional modulus.
    properties
        q {mustBeModulus(q)} = [];                   % modulus
        coefficients (1,:) {mustBeCoefficients} = 0; % highest power first
    end

    properties (Dependent)
        N (1,1)                                      % ring dimension
    end

    methods
        % constructor
        function c = RingElement(coefficients,q,options)
            % Instantiate a polynomial with the given coefficients. Optionally,
            % a modulus `q` can be specified.
            %
            % Other options may be passed via the struct `options` supporting
            % the following fields:
            %   N: integer
            %       Sets the number of coefficients explicitly.
            %   mode:
            %       How the unspecified higher-order coefficients are generated
            %       if `numel(coefficients)` is less than `options.N`. Either
            %       `"zeropad"` (insert zero coefficients) or `"repeat"` (insert
            %      copies of the highest-order coefficient of `coefficients`)
            %      are supported. The default is `"zeropad"`.
            arguments
                coefficients = [];
                q = [];
                options.N {mustBeInteger}
                options.mode {mustBeMember(options.mode,...
                   ["zeropad","repeat"])} = 'zeropad'
            end
            if nargin == 1
                c.coefficients = coefficients;
            elseif nargin == 2
                if ~isempty(q)
                    c.coefficients = Mod(coefficients,q);
                    c.q = q;
                else
                   c.coefficients = coefficients;
                end
            end
            if ~isempty(c.coefficients) && isfield(options,'N')
                if c.N < options.N
                    if strcmp(options.mode,'zeropad')
                        c.coefficients = [vpa(zeros(1,options.N-c.N)), ...
                            c.coefficients];
                    elseif strcmp(options.mode,'repeat')
                        c.coefficients = [repmat(c.coefficients(1),...
                            [1,options.N-c.N]), c.coefficients];
                    end
                end
            end
        end

        function N = get.N(c)
            % The number of coefficients
            N = length(c.coefficients);
        end

        % computations
        function c = uminus(c)
            % Negate all coefficients
            c.coefficients = -c.coefficients;
        end

        function c = plus(a,b)
            % Coefficient-wise addition
            [c,coefficients1,coefficients2] = handleinputs(a,b);
            c.coefficients = Mod(coefficients1+coefficients2,c.q);
        end

        function c = minus(a,b)
            % Coefficient-wise subtraction
            [c,coefficients1,coefficients2] = handleinputs(a,b);
            c.coefficients = Mod(coefficients1-coefficients2,c.q);
        end

        function c = times(a,b)
            % Coefficient-wise multiplication
            [c,coefficients1,coefficients2] = handleinputs(a,b);
            c.coefficients = Mod(coefficients1.*coefficients2,c.q);
        end

        function c = rdivide(a,b)
            % Coefficient-wise division
            [c,coefficients1,coefficients2] = handleinputs(a,b);
            c.coefficients = Mod(coefficients1./coefficients2,c.q);
        end

        function c = mtimes(a,b)
            % Polynomial multiplication. If one operand is scalar, the
            % multiplication is carried out as if the scalar was a constant
            % polynomial.
            [c,coefficients1,coefficients2] = handleinputs(a,b);
            if isscalar(coefficients1) || isscalar(coefficients2)
                c.coefficients = Mod(coefficients1*coefficients2,c.q);
            else
                c.coefficients = polymultRq(coefficients1,coefficients2,c.q);
            end
        end

        function c = mrdivide(a,b)
            % Coefficient-wise division. One of the operands must be scalar.
            [c,coefficients1,coefficients2] = handleinputs(a,b);
            if isscalar(coefficients1) || isscalar(coefficients2)
                c.coefficients = Mod(coefficients1./coefficients2,c.q);
            else
                error('polynomial inverse not defined')
            end
        end

        function c = power(a,b)
            % Coefficient-wise exponentiation of `a` to the power of `b`. Note
            % that `power(a, 0)` is defined as the polynomial with all
            % coefficients equal to one.
            arguments
                a {mustBeA(a,'RingElement')}
                b {mustBeScalarOrEmpty}
            end
            c = RingElement([],a.q);
            c.coefficients = Mod(a.coefficients.^b,a.q); % powermod would be faster
        end

        function [c,powa] = mpower(a,b)
            % (faster) Exponentiation of polynomial `a` to the power of `b`
            arguments
                a {mustBeA(a,'RingElement')}
                b {mustBeScalarOrEmpty}
            end
            bin_exp = floor(log2(b));
            powa = cell(bin_exp,1);
            powa{1} = a;
            for i = 1:bin_exp
                powa{i+1} = powa{i}*powa{i};
            end
            selection = de2bi(b-2^bin_exp,bin_exp);
            c = powa{end};
            for i = 1:length(selection)
                if selection(i) == 1
                    c = c*powa{i};
                end
            end
        end

        function c = innerproduct(a,b)
            % inner product of the coefficient vectors
            % if inputs are not of the same length, zeros are padded
            [c,coefficients1,coefficients2] = handleinputs(a,b);
            c.coefficients = Mod(coefficients1*coefficients2',c.q);
        end

        % rounding etc.
        function c = floor(c)
            % Coefficient-wise floor
            c.coefficients = floor(c.coefficients);
        end

        function c = ceil(c)
            % Coefficient-wise ceil
            c.coefficients = ceil(c.coefficients);
        end

        function c = round(c)
            % Coefficient-wise round
            c.coefficients = round(c.coefficients);
        end

        % modulus/scale
        function c = modreduce(c,qprime)
            % Reduce all coefficients modulo the new modulus `qprime`
            arguments
                c {mustBeA(c,'RingElement')}
                qprime {mustBeModulus(qprime)}
            end
            c.coefficients = Mod(c.coefficients,qprime);
            c.q = qprime;
        end

        function c = rescale(a,b)
            % Rescale `a` by `b`. This modifies the coefficients and the modulus.
           arguments
                a {mustBeA(a,'RingElement')}
                b {mustBePositive,mustBeScalarOrEmpty}
            end
            c = round(rdivide(a,b));
            c.q = round(c.q/b);
        end

        % evaluation
        function value = evalat(c,x)
            % Treat the ring element as a polynomial and evaluate it at the
            % given argument `x`.
            arguments
                c {mustBeA(c,'RingElement')}
                x {mustBeScalarOrEmpty}
            end
            value = vpa(0); % make type dependent
            for i = 0:c.N-1
                value = value + c.coefficients(end-i)*x^i;
            end
        end

        function c = evalpoly(a,poly)
            % Evaluate an expression that is polynomial in the ring element. The
            % multiplicative factors for each power of the ring element are
            % given by `poly` resulting in
            %   poly(1)*a.^(n-1) + poly(2)*a.^(n-2) + ... + poly(n)
            arguments
                a {mustBeA(a,'RingElement')}
                poly {mustBeCoefficients(poly)}
            end
            n = length(poly);
            c = RingElement(vpa(zeros(1,a.N)),a.q);
            for i = 1:n-1
                val = poly(i)*a.^(n-i);
                c = c + val;
            end
            % constant polynomial term is a special case because `a.^0` would
            % yield the polynomial with *all* coefficients equal to 1
            c = c + poly(n);
        end

        function value = norm(c,type)
            % Compute the norm of the vector of coefficients. For valid values
            % of `type`, refer to the built-in `norm` function.
            arguments
                c {mustBeA(c,'RingElement')}
                type
            end
            value = norm(c.coefficients,type);
        end

        % shape
        function c = rotate(c,positions)
            % Rotate the coefficients by the given number of positions. If
            % `positions` is positive, higher order monomial coefficients are
            % shifted towards lower orders, and vice versa if `positions` is
            % negative.
            arguments
                c {mustBeA(c,'RingElement')}
                positions {mustBeInteger}
            end
            c.coefficients = circshift(c.coefficients,positions);
            if positions < 0
                c.coefficients(end+positions+1:end) = -c.coefficients(end+positions+1:end);
            end
            if positions > 0
                c.coefficients(1:positions) = -c.coefficients(1:positions);
            end
        end

        function a = lt(a,b)
            % Short-hand for `rotate(a, -b)` with positive `b`.
            arguments
                a {mustBeA(a,'RingElement')}
                b {mustBeInteger,mustBePositive}
            end
            a = rotate(a,-b);
        end

        function a = gt(a,b)
            % Short-hand for `rotate(a, b)` with positive `b`.
             arguments
                a {mustBeA(a,'RingElement')}
                b {mustBeInteger,mustBePositive}
            end
            a = rotate(a,b);
        end

        function c = horzcat(a,b)
            % Stack the coefficients of `a` and `b` resulting in a new ring
            % element with `a.N` and `b.N` slots. `a`'s coefficients are used
            % for the higher-order monomials.
            [c,coefficients1,coefficients2] = handleinputs(a,b);
            c.coefficients = [coefficients1,coefficients2];
        end

        % display
        function disp(c,num_digits)
            if nargin == 1
                num_digits = 7; % default value
            end
            divisors = (floor(log10(abs(c.coefficients))) - num_digits + 1);
            divisors(divisors<1) = 0;
            first_digits = double(floor(c.coefficients ./ 10.^divisors));
            first_digits(divisors == 0) = c.coefficients(divisors == 0);
            divisors = double(divisors);

            if all(first_digits == 0)
                disp('0  (mod q) (mod X^N+1)');
                return
            else
                d = c.N - 1;
                s = cell(1,d);
                ind = 1;
                for i = 1:c.N
                    if first_digits(i) ~= 0
                        if ind ~= 1
                            if first_digits(i) > 0
                                s(ind) = {' + '};
                                ind = ind + 1;
                            else
                                s(ind) = {' - '};
                                first_digits(i) = -first_digits(i);
                                ind = ind + 1;
                            end
                        end
                        if first_digits(i) ~= 1 || d == 0
                            if first_digits(i) == -1
                                if i == c.N
                                    s(ind) = {'-1'};
                                else
                                    s(ind) = {'-'};
                                end
                                ind = ind + 1;
                            else
                                if divisors(i) ~= 0
                                    s(ind)={append(num2str(first_digits(i)),'...')};
                                else
                                    s(ind) = {num2str(first_digits(i))};
                                end
                                ind = ind + 1;
                                if d > 0
                                    s(ind) = {'*'};
                                    ind = ind + 1;
                                end
                            end
                        end
                        if d >= 2
                            s(ind) = {['X^' int2str(d)]};
                            ind = ind + 1;
                        elseif d == 1
                            s(ind) = {'X'};
                            ind = ind + 1;
                        end
                    end
                    d = d - 1;
                end
            end
            disp([[s{:}],'  (mod q) (mod X^N+1)']);
        end

    end

    methods (Hidden)
        % input handling for computation
        function [c,coefficients1,coefficients2] = handleinputs(a,b)

            if isa(a,'RingElement') && isa(b,'RingElement')
                if a.q < b.q
                    coefficients1 = a.coefficients;
                    coefficients2 = Mod(b.coefficients,a.q);
                    warning('modulus of second input is reduced')
                    c = RingElement([],a.q);
                elseif a.q > b.q
                    coefficients1 = Mod(a.coefficients,b.q);
                    coefficients2 = b.coefficients;
                    warning('modulus of first input is reduced')
                    c = RingElement([],b.q);
                else % equal moduli
                    coefficients1 = a.coefficients;
                    coefficients2 = b.coefficients;
                    c = RingElement([],a.q);
                end
                if a.N < b.N
                    coefficients1 = [vpa(zeros(1,b.N-a.N)) coefficients1];
                elseif b.N < a.N
                    coefficients2 = [vpa(zeros(1,a.N-b.N)) coefficients2];
                end

            elseif isa(a,'RingElement') && ~isa(b,'RingElement')
                if ~isnumeric(b) && ~isa(b,'sym')
                    error('Second argument must be a RingElement or numeric')
                end
                if isscalar(b)
                    coefficients2 = b;
                elseif isvector(b)
                    if length(b) < a.N
                        coefficients2 = [vpa(zeros(1,a.N-length(b))) b];
                        coefficients2 = reshape(coefficients2,[1 a.N]);
                    elseif length(b) > a.N
                        error('vector dimension too big')
                    else
                        coefficients2 = reshape(b,[1 a.N]);
                    end
                end
                coefficients1 = a.coefficients;
                c = RingElement([],a.q);

            elseif isa(b,'RingElement') && ~isa(a,'RingElement')
                if ~isnumeric(a) && ~isa(a,'sym')
                    error('Second argument must be a RingElement or numeric')
                end
                if isscalar(a)
                    coefficients1 = a;
                elseif isvector(a)
                    if length(a) < b.N
                        coefficients1 = [vpa(zeros(1,b.N-length(a))) a];
                        coefficients1 = reshape(coefficients1,[1 b.N]);
                    elseif length(a) > b.N
                        error('vector dimension too big')
                    else
                        coefficients1 = reshape(a,[1 b.N]);
                    end
                end
                coefficients2 = b.coefficients;
                c = RingElement([],b.q);
            else
                error('this case should be impossible')
            end
        end
    end

end