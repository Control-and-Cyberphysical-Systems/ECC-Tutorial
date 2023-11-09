function [zvec,z] = Dcd(m,Encoder)
    arguments
        m
        Encoder {mustBeEncoder(Encoder)} = [];
    end

    is_scalar = 0;
    if ~iscell(m) % scalar case
        is_scalar = 1;
        mustBeA(m,'RingElement')
        if isempty(m.s) || isempty(m.type)
            error('decoding failed due to unspecified encoding parameters.')
        end
    end

    if is_scalar
        if m.type == "packed"
            if isempty(Encoder)
                Encoder = CoefficientEncoder(m.N/2,m.type);
            end
            z = Encoder.decode(m,m.s);
            zvec = z.coefficients;
        end
        if m.type == "scalar"
            if isempty(Encoder)
                Encoder = CoefficientEncoder(m.N,m.type);
            end
            z = Encoder.decode(m,m.s);
            zvec = z.coefficients;
        end
    end

    if ~is_scalar % matrix case
        if m{1,1}.type == "packed" && isempty(Encoder)
            Encoder = CoefficientEncoder(m{1,1}.N/2,m{1,1}.type);
        elseif m{1,1}.type == "scalar" && isempty(Encoder)
            Encoder = CoefficientEncoder(m{1,1}.N,m{1,1}.type);
        end
        zvec = [];
        for i = 1:size(m,1)
            for j = 1:size(m,2)
                mustBeA(m{i,j},'RingElement')
                z{i,j} = Encoder.decode(m{i,j},m{i,j}.s);
                zvec = [zvec; z{i,j}.coefficients];
            end
        end
    end
end