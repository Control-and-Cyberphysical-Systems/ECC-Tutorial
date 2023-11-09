classdef Setup
    % Bundles the configuration
    properties
        N {mustBePositive,mustBeInteger}           % key dimension
        q {mustBeModulus(q)} = [];                 % modulus
        r {mustBePositive}                         % error bound
    end

    methods
        % constructor
        function setup = Setup(N,q,r)
            arguments
                N
                q
                r
            end
            setup.N = N;
            setup.q = q;
            setup.r = r;
        end
    end

end