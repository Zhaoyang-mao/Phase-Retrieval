function x0 = LS_init(Yhat, W, L)
% This function implements the initialization method for phase retrieval in STFT.
%�ú���ʵ������λ�ָ���STFT��ʼ���ķ���
% The initialization is based on the principal eigenvector (with appropriate
% normalization) of a matrix, constructed as the solution of a least-squares
% problem.
%
% Inputs:
%  Yhat: FFT of the squared absolute values of the STFT of a signal (noisy)
%Yhat���ź�STFT����ֵƽ����FFT
%  W   : window support length
%W:�������ĳ���
%  L   : separation in time between adjacent sections
%L�����������ʱ����
% Outputs:
%  x0  : initial estimation of the signal
%x0:�źŵĳ�ʼ����
% Written by Tamir Bendory; last update 28/07/2016.


	%% ------------------ PRELIMINARY CALCULATIONS ��������

	% signal's length
	N = size(Yhat,2); vec=(0:N-1)';     %�����źų��ȣ�

	% STFT rectangular window           %���ɾ��δ�
	g = zeros(N,1); g(1:W)=1;

	%%------ CONSTRUCTING THE APPROXIMATION MATRIX-----------------
            %������ƾ���
	x_est = zeros(N);

	for ii = -(W-1):W-1
		
		% interpolating the missing entries of y with a cubic interpolation
		y = interp1(1:L:N, Yhat(:, mod(ii, N)+1), 1:N, 'pchip').';
		
		% recall that y=Gx=F*\SigmaFx. So, we aim at solving 
		% min_x (\|Fy - \SigmaFx\|)^2 for each ii. This is performed by the fft.
		
		S = diag(fft(g.*g(mod(vec-ii,N)+1)));
		xhat_est = S\fft(y);
		x_est(:, mod(ii, N)+1) = ifft(xhat_est);
		
	end

	%% ----------- SIGNAL ESTIMATION �źŹ���-------------------------

	% reorder the solution as Xd=x*x';
	xd = convert2diag(x_est);

	%extracting the principal eigenvector of the approximation matrix
    %��ȡ�������������
	[eigvec, eigval] = eigs(xd,1);

	% normalization  ��׼��
	d = diag(xd);
	x0 = sqrt(sum(  d(d>0)  ))*eigvec;

end
