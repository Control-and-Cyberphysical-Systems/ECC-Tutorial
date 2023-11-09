function sk = KeyGen(setup)
% Generates the secret key
    arguments
        setup {mustBeA(setup,'Setup')}
    end
    % NOTE: The randomness and emulated distributions do not satisfy
    %       cryptographic requirements (see comments below)!
    sk = Mod(round(rand([setup.N, 1])*setup.q),setup.q);
end