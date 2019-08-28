function x = my_stft_i(X, L, g)
% Pseudo-inverse of the STFT transform
% Inputs:
%   X: STFT matrix
%   L: separation in time between adjacent sections
%   g: the analysis STFT window
% Outputs:
%   x: signal whose STFT is closest to X in least squares sense
%
% Authors: Tamir Bendory and Nicolas Boumal, July 2017
%
% See also: my_stft

    N = length(g);
    gmn = zeros(ceil(N/L),N);
    vec = 0:(N-1);
    vecL = 0:L:(N-1);

    for ii = 0:ceil(N/L)-1
        gmn(ii+1,:) = g(mod(-vec+vecL(ii+1),N)+1);
    end
    Gabs = sum(abs(gmn).^2, 1);
    %X_szie = size(X);
    iFX = (ifft(X.')).';%iFX size Ϊ15*29
    x = (sum(iFX.*conj(gmn), 1)./Gabs).';

end
