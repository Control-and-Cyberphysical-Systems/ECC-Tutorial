function [z,Encoder] = Ecd(setup,x,options)
% A wrapper for CoefficientEncoder
    arguments
        setup {mustBeA(setup,'Setup')}
        x {mustBeCoefficients(x)}       
        options.type {mustBeMember(options.type,...
                   ["scalar","packed"])} = 'packed'
        options.mode {mustBeMember(options.mode,...
                   ["single","zeropad","repeat"])} = 'single'
    end

    if options.type == "packed" && size(x,2)*2>setup.N
        error('packed encoding can only encode up to N/2 coefficients')
    end
    if options.type == "scalar" && size(x,2)>setup.N
        error('scalar encoding can only encode up to N coefficients')
    end

    Encoder = CoefficientEncoder(setup.N/2,options.type);
    
    % rows of x get encoded in polynomial coefficients
    if options.type == "packed"
        for i = 1:size(x,1)
                p = RingElement(vpa(x(i,:)),setup.q0);
                z{i} = Encoder.encode(p,setup.s);
                z{i}.s = setup.s;
                z{i}.type = options.type;
                z{i}.mode = options.mode;
        end
        if size(x,1) == 1
            z = z{1};
        end

    end

    if options.type == "scalar"
        if options.mode == "single"
            for i = 1:size(x,1)
                for j = 1:size(x,2)
                    p = RingElement(vpa(x(i,j)),setup.q0,...
                    'N',setup.N,'mode','zeropad');
                    z{i,j} = Encoder.encode(p,setup.s);
                    z{i,j}.s = setup.s;
                    z{i,j}.type = options.type;
                    z{i,j}.mode = options.mode;
                end
            end
            if size(x,1)== 1 && size(x,2) == 1
                z = z{1,1};
            end
        elseif options.mode == "zeropad" || options.mode == "repeat"
            for i = 1:size(x,1)
                p = RingElement(vpa(x(i,:)),setup.q0,...
                    'N',setup.N,'mode',options.mode);
                z{i} = Encoder.encode(p,setup.s);
                z{i}.s = setup.s;
                z{i}.type = options.type;
                z{i}.mode = options.mode;
            end
            if size(x,1) == 1
                z = z{1};
            end
        end
    end
end