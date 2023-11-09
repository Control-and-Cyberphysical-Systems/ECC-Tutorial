function [a,b] = adaptModScale(a,b,mode,scaling)
% adapts the modulus or the scaling such that a,b are compatible
% modreduce -> reduces modulus keeps scaling
% rescale   -> reduces modulus and scaling

if mode == "auto"
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
elseif mode == "manual" % adapt manually
    if a.q ~= b.q
        error("incompatible moduli." + ...
            " use modulus reduce or set adapt to 'auto'")
    elseif a.s ~= b.s
        error("incompatible scaling." + ...
            "use rescale or set adapt to 'auto'")
    end
end

end