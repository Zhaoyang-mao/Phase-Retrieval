clear
%%======================================================================%%
%本历程是复现Yonina C Eldar的文章 Sparse Phase Retrieval from Shot-Time
%Fourier Measurements 中的Griffin-Lim algorithm(GLA)算法
%%======================================================================%%
N = 29;                         %设置信号的长度%这里的N与W应满足N>=2W-1并且
                                %N和W-1为互为素数
W = 12 ;                         %设置窗的大小，该脚本为一个矩形窗
L = 2;                          %设置相邻窗之间的时间间隔
%%
%产生真值信号
% x = randn(N,1);                 %真值信号为N*1维的矩阵
%%
%产生离散chirp信号
k=-14:14;
N = length(k);
x = real((exp(1i*0.1141.*k.^(2)))');
%%
%产生高斯窗
g = zeros(N,1);                 %产生窗函数，先生成一个0矩阵
% guss_x = 0:W-1;
% guss_sigma = 2;                 %调节高斯窗的sigma
% guss_w = normpdf(guss_x,W/2+1,guss_sigma);      %生成离散的高斯窗信号
% %g(1:W) = 1;                     %让0矩阵的前W项为1，也就是窗的大小
% g(1:W)=guss_w;
g1=gausswin(W);          %高斯窗
% g1 = hanning(W);       %海明窗
% g1 = hann(W);             %汉宁窗
% g1 = blackman(W);          %布莱克曼窗
g(1:W) = g1;
%%
%产生矩形窗
% g = zeros(N,1);                 %产生窗函数，先生成一个0矩阵
% g(1:W) = 1;                     %使得0矩阵的前W项为1，也就是窗的大小
%%
X = my_stft(x,L,g);             %进行STFT变换,也就是文章中的式1
Yc = abs(X).^2;                 %取出振幅的平方，也就是文中算法中的第一步

%%===============噪音的加入============================================
snr = 20;                        %信噪比，越高噪声越小，越低噪音越大
noise = randn(size(Yc));        %生成一个与Yc大小相同的随机矩阵
sigma = norm(Yc,'fro')*10^(-snr/20)/norm(noise,'fro');
Y = Yc+sigma*noise;             %将噪声融入到将要输入的信号中

T = tic();                      %启动秒表计时器
Yhat = fft(Y')'/N;       
x0 = LS_init(Yhat,W,L);         %将初始处理的信号x0进行初始化

x_hat_gla = PR_STFT_GLA(Y,N,L,g,x0);
x_hat_gla = real(x_hat_gla*sign(x_hat_gla'*x));

%=======================================================================%
%绘图
%=======================================================================%

my_plot(x,x_hat_gla);
