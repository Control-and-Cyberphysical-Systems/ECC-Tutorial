function ab = polymultRq(a,b,q,options)
%% Polynomial multiplication in Rq=Zq[x]/(X^N+1), where X^N+1 is irreducible
% Input: a - vector of polynomial coefficients from a=sum_{i=0}^{N-1} a[i] * X^i
%        b - vector of polynomial coefficients from b=sum_{i=0}^{N-1} b[i] * X^i
%        q - ciphertext modulus
%
% Output: a(X)*b(x) (mod q) (mod X^N+1) - polynomial coefficients in Rq
%
% Methods: 1) polynomial multiplication + polynomial division ('naive')
%          2) Negative wrapped convolution ('conv')
arguments
    a {mustBeCoefficients(a)}
    b {mustBeCoefficients(b)}
    q {mustBeModulus(q)}
    options.mode {mustBeMember(options.mode,["naive","conv","ntt"])} = 'conv';
end

switch options.mode
    case 'naive' % convolution and division
        N = length(a);
        ab = polymult(a,b,q);
        R = vpa([1 zeros(1,N-1) 1]); % [X^N, 0,...,0,1]
        [~,ab] = polydiv(ab,R,q);

    case 'conv' % O(n^2) complexity (no division)
        ab = negativewrappedconvolution(a,b,q);

    % case 'ntt'
    % not yet implemented
end

end