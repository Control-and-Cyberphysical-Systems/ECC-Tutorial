function w = polymult(a,b,q)
%% naive implementation of polynomial multiplication
% using convolution w(k)=sum_j u(j)*v(k-j+1) in matrix form

m=length(a);
n=length(b); 
lenw=m+n-1;

U=[zeros(1,m-1),a,zeros(1,m-1)];

A=vpa(zeros([n,lenw]));
for i=1:lenw
   A(:,i)=fliplr(U(i:i+n-1));
end

w=Mod(b*A,q);

end