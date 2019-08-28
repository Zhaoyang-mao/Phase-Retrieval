function [S,F,T,P] = STFT_M(x,varargin)
%function [S,F,T,P] = STFT_M(x,window,noverlap,nfft,fs)为函数原型
%可变长度输入varargin从该点开始收集输入的参数
%=========================================================================%
%功能：使用短时傅里叶变换得到信号的频谱图
%说明;在使用时无参数输出，会自动绘制频谱图；有参数输出，则会返回输入信号的短时傅里
%变换，当然也可以从函数的返回值S F T P绘制频谱图
%%%%%%%%%%%%参数描述
%x          输入信号的向量，默认情况下，即x之后的参数没有输入时，x将被分为8段做变换处
%           理，如果x不能被平分成8段，则会做截断处理
%window     窗函数，默认情况下为nfft长度的海明窗Hamming。如果window为一个证书，x
%           将被分成x段，每段使用海明窗加窗。如果window是一个向量，x将被分为
%           length(window)段，每一段使用window向量制定的窗函数加窗
%noverlap   每一段的重叠样本数，默认情况是在各段之间产生50%的重叠。它必须为一个
%           小于window或者length(widow)的整数。其意思为两个相邻窗不是尾接头的，
%           而是两窗之间是有交集的，有重叠的部分
%nfft       做fft变换的长度，默认情况下维256和大于每段长度最小2次幂之间的最大值
%           另外，此参数除使用一个常量外，还可以指定一个频率向量F。它需要为标量
%fs         采样频率，默认值归一化频率，如果制定为[],默认为1Hz
%%%%%%%%%%%%返回值描述
%S          输入信号x的短时傅里变换。它的每一列包含一个短期局部时间的频率成分估计
%           时间沿列增加，频率沿行增加。如果x是长度为Nx的复信号，则S为nfft行k列
%           的复矩阵，其中k取决于window，如果window为一个标量，则
%           k = fix((Nx-noverlap)/(window-noverlap))；如果window为向量，则
%           k = fix((Nx-noverlap)/(length(window)-noverlap))。对于实信号x，
%           如果nfft为偶数，则S的行数为(nfft/2+1)，如果nfft为奇数，则行数为
%           (nfft+1)/2，列数同上
%F          在输入变量中使用F频率向量，函数会使用Goertzel方法计算在F指定的频率
%           处计算频谱图。指定的频率被四舍五入到与信号分辨率相关的最近的DFT
%           容器(bin)中。而在其他的使用nfft语法中，短时傅里叶变换方法将被使用。
%           对于返回值中的F向量，为四舍五入的频率，其长度等于S的行数。
%T          频谱图计算的时刻点，其长度等于上面定义的k值，值为所分各段的中点
%P          能量谱密度PSD(Power Spectral Density)，对于实信号，P是各段PSD的单
%           边周期估计；对于复信号，当指定F频率向量时，P为双边PSD。P矩阵的元素
%           计算公式如下P(I,j)=k|S(I,j)|2，其中的的k是实值标量，定义如下对于单
%           边PSD，计算公式如下，其中w(n)表示窗函数，Fs为采样频率，在0频率和奈
%           奎斯特频率处，分子上的因子2改为1；
%关于画图的问题，如果采样频率没有指定，fs由2*pi来代替，当调用函数时，没有输出参数
%将会自动绘制各段的PSD估计，绘制命令如下
%surf(T,F,10*log10(abs(p)));
%axis tight;
%view(0,90);
%注意事项：nfft越大，频域的分辨率就越高（分辨率=fs/nfft）,但是离瞬时频率就越远
%novelap影响时间轴的分辨率，越接近bfft，分辨率就越高，相应的冗余就越多，计算量越大
%=========================================================================%
[S,F,T,P] = spectrogram(x,varargin);