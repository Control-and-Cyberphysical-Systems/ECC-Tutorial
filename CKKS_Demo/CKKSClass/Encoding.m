classdef (InferiorClasses = {?sym,?RingElement}) Encoding < RingElement
    % A RingElement extended by scaling and type properties.
    % See the base class documentation.
    properties
        s    {mustBeScalingFactor(s)} = [];        % scaling factor
        adapt {mustBeAdaptionMode(adapt)} = 'auto' % modulo and scaling adaptation
    end
    properties (SetAccess = ?Encoder)
        type {mustBeEncodingType(type)}            % encoding type
    end

    methods
        % constructor
        function m = Encoding(coefficients,q,s,type,adapt)
            % This constructor is only used internally. See
            % `CoefficientEncoder`, for instance.
            if nargin >= 1
                m.coefficients = coefficients;
            end
            if nargin >= 2
                m.q = q;
                m.coefficients = Mod(m.coefficients,q);
            end
            if nargin >= 3
                m.s = s;
            end
            if nargin >= 4
                m.type = type;
            end
            if nargin == 5
                m.adapt = adapt;
            end
        end

        % computation (subclass outputs RingElements)
        function c = plus(a,b)
            % See documentation of `RingElement.plus`
            [a,b,cs,ctype,cadapt] = handleinputs(a,b,'constant');
            c = plus@RingElement(a,b);
            c = Encoding(c.coefficients,c.q,cs,ctype,cadapt);
        end

        function c = minus(a,b)
            % See documentation of `RingElement.minus`
            [a,b,cs,ctype,cadapt] = handleinputs(a,b,'constant');
            c = minus@RingElement(a,b);
            c = Encoding(c.coefficients,c.q,cs,ctype,cadapt);
        end

        function c = times(a,b)
            % See documentation of `RingElement.times`
            [a,b,cs,ctype,cadapt] = handleinputs(a,b,'increase');
            c = times@RingElement(a,b);
            c = Encoding(c.coefficients,c.q,cs,ctype,cadapt);
        end

        function c = rdivide(a,b)
            % See documentation of `RingElement.rdivide`
           [a,b,cs,ctype,cadapt] = handleinputs(a,b,'increase');
            c = rdivide@RingElement(a,b);
            c = Encoding(c.coefficients,c.q,cs,ctype,cadapt);
        end

        function c = mtimes(a,b)
            % See documentation of `RingElement.mtimes`
            [a,b,cs,ctype,cadapt] = handleinputs(a,b,'increase');
            c = mtimes@RingElement(a,b);
            c = Encoding(c.coefficients,c.q,cs,ctype,cadapt);
        end

        function c = mrdivide(a,b)
            % See documentation of `RingElement.mrdivide`
            [a,b,cs,ctype,cadapt] = handleinputs(a,b,'increase');
            c = mrdivide@RingElement(a,b);
            c = Encoding(c.coefficients,c.q,cs,ctype,cadapt);
        end

        function c = power(a,b)
            % See documentation of `RingElement.power`
            cs = a.s;
            ctype = a.type;
            [a,b,~,~,cadapt] = handleinputs(a,b,'increase');
            c = power@RingElement(a,b);
            c = Encoding(c.coefficients,c.q,cs^b,ctype,cadapt);
        end

        function c = mpower(a,b)
            % See documentation of `RingElement.mpower`
            cs = a.s;
            ctype = a.type;
            [a,b,~,~,cadapt] = handleinputs(a,b,'increase');
            c = mpower@RingElement(a,b);
            c = Encoding(c.coefficients,c.q,cs^b,ctype,cadapt);
        end

        function c = horzcat(a,b)
            % See documentation of `RingElement.horzcat`
            [a,b,cs,ctype,cadapt] = handleinputs(a,b,'constant');
            c = horzcat@RingElement(a,b);
            c = Encoding(c.coefficients,c.q,cs,ctype,cadapt);
        end

        function c = innerproduct(a,b)
            % See documentation of `RingElement.innerproduct`
            [a,b,cs,ctype,cadapt] = handleinputs(a,b,'increase');
            c = innerproduct@RingElement(a,b);
            c = Encoding(c.coefficients,c.q,cs,ctype,cadapt);
        end

        % modulus/scale
        function c = rescale(a,b)
            [a,b,cs,ctype,cadapt] = handleinputs(a,b,'reduce');
            c = rescale@RingElement(a,b);
            c = Encoding(c.coefficients,c.q,cs,ctype,cadapt);
        end

        % RS variants (keeping the scaling constant)
            % Variant of `to,es` with subsequent rescaling
        function c = RStimes(a,b)
            c = times(a,b);
            c = rescale(c,c.s/a.s);
        end

        function c = RSrdivide(a,b)
            % Variant of `rdivide` with subsequent rescaling
            c = rdivide(a,b);
            c = rescale(c,c.s/a.s);
        end

        function c = RSmult(a,b)
            % Variant of `mtimes` with subsequent rescaling
            c = mtimes(a,b);
            c = rescale(c,c.s/a.s);
        end

        function c = RSmrdivide(a,b)
            % Variant of `mrdivide` with subsequent rescaling
            c = mrdivide(a,b);
            c = rescale(c,c.s/a.s);
        end

        function c = RSpow(a,b)
            % Variant of `power` with subsequent rescaling
            c = power(a,b);
            c = rescale(c,c.s/a.s);
        end

        function c = RSmpower(a,b)
            % Variant of `mpower` with subsequent rescaling
            c = mpower(a,b);
            c = rescale(c,c.s/a.s);
        end

        function c = RSinnerproduct(a,b)
            c = innerproduct(a,b);
            c = rescale(c,c.s/a.s);
        end

        % display
        function disp(c)
            disp@RingElement(c)

            if isempty(c.s) && isempty(c.type)
                disp('unspecified encoding with unspecified scaling')
            elseif ~isempty(c.s) && isempty(c.type)
                disp(['unspecified encoding with scaling = ',num2str(double(c.s))])
            elseif isempty(c.s) && ~isempty(c.type)
                disp([c.type,'-encoding with unspecified scaling'])
            elseif ~isempty(c.s) && ~isempty(c.type)
                disp([c.type,'-encoding with s = ',num2str(double(c.s))])
            end
        end

    end

    methods (Hidden)
        function [a,b,cs,ctype,cmode] = handleinputs(a,b,scaling)
            if isa(a,'Encoding') && isa(b,'Encoding')
                if ~strcmp(a.type,b.type)
                    error('both inputs must share the same encoding type')
                end
                [a,b,cmode] = adaptModScale(a,b,scaling);
                % compute the scaling
                if strcmp(scaling,'increase')
                    cs = a.s*b.s;
                elseif strcmp(scaling,'constant')
                    cs = a.s;
                end
                ctype = a.type;
                a = RingElement(a.coefficients,a.q);
                b = RingElement(b.coefficients,b.q);
            elseif isa(a,'Encoding') && ~isa(b,'Encoding')
                cmode = a.adapt;
                if strcmp(scaling,'reduce')
                    cs = round(a.s/b);
                else
                    cs = a.s;
                end
                ctype = a.type;
                a = RingElement(a.coefficients,a.q);
            elseif ~isa(a,'Encoding') && isa(b,'Encoding')
                cs = b.s;
                ctype = b.type;
                cmode = b.adapt;
                b = RingElement(b.coefficients,b.q);
            end
        end
    end
end