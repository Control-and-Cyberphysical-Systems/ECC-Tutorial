function [y,r] = polydiv(b,a,q)
%% naive implementation of polynomial division
% computes y(x),r(x) given b(x),a(x) s.t. b(x) = a(x)*y(x) + r(x)

n = length(b)-1;   
m = length(a)-1;
a = vpa([a,zeros(1,n-m)]);
y = vpa(0);

for  k = 1:n+1
    if k < n-m+2
        y(k) = (b(k)-[y(1:k-1)]*[a(k:-1:2)].') /a(1);
    else
        r(k-(n-m+1)) = b(k)-[y(1:n-m+1)]*[a(k:-1:k-n+m)].';
    end
end

if m-length(r)>0
    r=[zeros(1,m+1-length(r)),r]; % 0-padding
end

%y = Mod(y,q);
r = Mod(r,q);

end