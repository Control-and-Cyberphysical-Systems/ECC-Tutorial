classdef CoefficientEncoder < dynamicprops
    properties (SetAccess = immutable)
        M    {mustBePositive,mustBeInteger} = []; % dimension
        type {mustBeEncodingType(type)}           % encoding type
    end

    methods
        % initialize encoder
        function obj = CoefficientEncoder(N,type)
            % Instantiate an encoder.
            %
            % Inputs:
            %   N: integer
            %       The number of coefficients to encode.
            %   type:
            %       Either `'scalar'` or `'packed'`. The former is the identity.
            %       The latter involves the DFT-like transformation
            %       characteristic to CKKS.
            arguments
                N {mustBePositive,mustBeInteger}
                type {mustBeEncodingType(type)}
            end

            obj.type = type;
            switch type
                case 'scalar'
                    obj.M = N;   % all slots can be used
                case 'packed'
                    obj.M = 2*N; % half of the slots can be used
                    [U,CRT,CRTinv] = initializeEmbedding(obj,obj.M);
                    tmp = addprop(obj,'U');
                    obj.U = U;
                    tmp.SetAccess = 'immutable';
                    tmp = addprop(obj,'CRT');
                    obj.CRT = CRT;
                    tmp.SetAccess = 'immutable';
                    tmp = addprop(obj,'CRTinv');
                    obj.CRTinv = CRTinv;
                    tmp.SetAccess = 'immutable';
            end
        end

        % precompute canonical embedding
        function [U,CRT,CRTinv] = initializeEmbedding(~,M)
            if M < 2 || abs(floor(M/2)-M/2) > eps
                error('ring dimension must be a positive power of 2')
            end
            zeta = vpa(exp(-1i*pi/M));
            zeta_vec = zeta.^(1+4*(0:M-1)); % roots of unity
            U = vander(zeta_vec);
            U = U(1:M/2,:);
            CRT = [U; conj(U)];
            CRTinv = 1/M*conj(CRT.');       % embedding of C^(N/2) into polynomial ring
        end

        function [m,m_vec,z_ecd] = encode(obj,z,s)
            % Encode a RingElement object `z`
            arguments
                obj
                z {mustBeA(z,'RingElement')}
                s {mustBeScalingFactor}
            end
            if isempty(z.q)
                error('modulus must be specified');
            elseif isempty(s)
                error('scaling factor must be specified');
            end
            coefficients = z.coefficients;
            q = z.q;
            N = z.N;

            switch obj.type
                case 'scalar'
                    m_vec = round(s*coefficients);
                case 'packed'
                    if obj.M < 2*N
                        error('CKKS encoding supports N/2 slots at most')
                    elseif obj.M > 2*N
                        if isa(coefficients,'sym')
                            padding = vpa(zeros(1,obj.M/2-N));
                        else
                            padding = zeros(1,obj.M/2-N);
                        end
                        coefficients = [padding, coefficients];

                    end
                    lift = [coefficients(:); conj(coefficients(:))];
                    z_ecd = obj.CRTinv*lift;
                    m_vec = round(s*obj.CRTinv*lift);
            end
            m = RingElement(m_vec,q);
        end

        function [polynomial,coefficients,z_dcd] = decode(obj,m,s)
            % Decoding the Encoding object `m` to a corresponding RingElement object.
            arguments
                obj
                m
                s {mustBeScalingFactor(s)}
            end
            switch obj.type
                case 'scalar'
                    if obj.M ~= m.N
                        error('encoder and polynomial do not share dimension')
                    end
                    coefficients = m.coefficients/s;
                case 'packed'
                    z_dcd = obj.U*m.coefficients(:);
                    coefficients = real(obj.U*m.coefficients(:)/s)';
            end
           polynomial = RingElement(coefficients,m.q);
        end
    end

end