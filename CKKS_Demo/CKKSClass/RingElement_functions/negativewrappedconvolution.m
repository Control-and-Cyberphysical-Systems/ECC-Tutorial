function ab = negativewrappedconvolution(a,b,q)
%% negacyclic convolution with MSB first
% a=[1 2 3 4];
% b=[0 0 0 5];
% A=[a(4)  a(3)   a(2)  a(1);
%   -a(1)  a(4)   a(3)  a(2);
%   -a(2) -a(1)   a(4)  a(3);
%   -a(3) -a(2)  -a(1)  a(4)];

N = length(a);
a = fliplr(a);
A = vpa(zeros(N,N));
for i = 1:N
    if i == 1
        A(i,:) = a;
    else
        a = circshift(a,1);
        a = [-a(1) a(2:end)];
        A(i,:) = a;
    end
end
ab = Mod(b*A',q);

end