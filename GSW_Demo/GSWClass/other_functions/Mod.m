function y = Mod(x,q)
    y = mod(x,q); 
    y = y - double(y>q/2)*q;
end