function [a,b,adapt_mode] = adaptModScale(a,b,scaling)
% adapts the modulus or the scaling such that a,b are compatible
% modreduce -> reduces modulus keeps scaling
% rescale   -> reduces modulus and scaling

mode = (strcmp(a.adapt,'auto') || strcmp(b.adapt,'auto'));

if mode == 1 % adapt automatically
    adapt_mode = 'auto';
    if strcmp(scaling,'increase')
        if a.q < b.q
            b = modreduce(b,a.q);
        elseif b.q < a.q
            a = modreduce(a,b.q);
        end
    elseif strcmp(scaling,'constant')
        if a.s < b.s
            b = rescale(b,b.s/a.s);
        elseif b.s < a.s
            a = rescale(a,a.s/b.s);
        end
        if a.q < b.q
            b = modreduce(b,a.q);
        elseif b.q < a.q
            a = modreduce(a,b.q);
        end
    end
elseif mode == 0 % adapt manually
    adapt_mode = 'manual';
    if a.q ~= b.q
        error("incompatible moduli." + ...
            " use modulus reduce or set adapt to 'auto'")
    elseif a.s ~= b.s
        error("incompatible scaling." + ...
            "use rescale or set adapt to 'auto'")
    end
end

end