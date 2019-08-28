function y = my_stft(x, L, g)
% Forward short-time Fourier transform (STFT)
% Inputs:
%   x: signal in the time domain
%   L: separation in time between adjacent sections
%   g: the analysis STFT window
% Outputs:
%   X: STFT of the signal (matrix)
%
% Authors: Tamir Bendory and Nicolas Boumal, July 2017
%
% See also: my_stft_i

	%% ---------------GENERATING PRELIMINARY VARIABLS----------------------

	N = length(x);
	vec = (0:(N-1))';
	m = (0:L:(N-1))';
	y = zeros(length(m), N);


	%% ----------- COMPUTES THE STFT

	for ii = 1:length(m)
		
		xg = x .* g(mod(m(ii) - vec, N) + 1);
		y(ii, :) = fft(xg);
		
	end

end
