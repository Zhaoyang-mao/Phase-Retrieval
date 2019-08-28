function xg = PR_STFT_GLA(Y, N, L, g, x0)
% GLA algorithm for phase retrieval from STFT measurements
%
% inputs:
% Y  - STFT magnitudes
% N - signal's length
% L - spacing between adjacent windows.
% g - stft window
% x0 - initial estimation
%
% output: estimated signal
%
% Authors: Tamir Bendory and Nicolas Boumal, July 2017


	gmn = zeros(ceil(N/L), N);
	vec = 0:(N-1);
	vecL = 0:L:(N-1);

    for ii = 0:ceil(N/L)-1
		gmn(ii+1,:) = g(mod(-vec+vecL(ii+1), N)+1);
    end

    
	err_iter = inf;
	iter = 1;

	xg = x0;
	max_iter = 5000;

    while err_iter > 1e-10 && iter < max_iter
		
		yg = my_stft(xg, L, g); % computing the STFT of xest
		b = sign(yg).*sqrt(Y); % keeping the phase  
		
		xg_new = my_stft_i(b, L, g); % inverting
		err_iter = norm(xg_new - xg, 2) / norm(xg, 2);
		xg = real(xg_new);
		iter = iter+1;
		
    end

end
