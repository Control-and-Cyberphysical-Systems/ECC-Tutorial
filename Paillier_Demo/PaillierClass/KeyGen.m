function [pk,sk] = KeyGen(lambda)
    arguments
        lambda {mustBePositive,mustBeInteger} = 1028
    end
    p1_bits = [1,datasample([0,1],lambda-1)];
    p2_bits = [1,datasample([0,1],lambda-1)];
    p1 = vpa(0);
    p2 = vpa(0);
    two = vpa(2);
    for i = 0:lambda-1
        p1 = p1 + p1_bits(end-i)*two^i;
        p2 = p2 + p2_bits(end-i)*two^i;
    end
    
    p1 = prevprime(p1);
    p2 = prevprime(p2);
    
    pk = p1*p2;         % should have 2048 binary places for security
    sk = (p1-1)*(p2-1); % gcd(pk,sk) = 1 ensured
end