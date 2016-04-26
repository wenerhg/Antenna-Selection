%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        K-user IA simulation       %
%           Xinyu Zhang             %
%            March 2014             %
%  Dalian University of Technology  %
%             ver 1.0               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
% close all

Mr = 2;    %接收天线数
Mt = 2;    %发射天线数
d = 1;    %数据流
K = 3;    %用户数
loop =  10; %循环次数
iter_leakmin = 50; %迭代步数

SNR = [0:5:100]; %发射功率范围
tic
 
for k = 1:loop     %%循环loop次取平均
    
    H = sqrt(1/2).*((normrnd(0,1,[Mr,Mt,K,K])+sqrt(-1).*normrnd(0,1,[Mr,Mt,K,K])));  %%产生独立同分布的循环对称复高斯分布的随机变量
    U = normalize(randn(Mr,d,K)); % 初始化干扰抑制矩阵
    V = normalize(randn(Mt,d,K)); % 初始化波束成形矩阵
    for snr = 1:length(SNR)
        
    xigma = 10^(SNR(snr)/10)/(1*d);  
    [V,U] = maxSINR_K_user(H, V, xigma, iter_leakmin); %迭代计算U，V
    R(snr,k) = rate_K_user(U, H, V, xigma);%计算总速率（bits/s/Hz)
    end
  k
end
display(['Total Time :']);
toc
R_1 = mean(R,2);
handles=plot(SNR,R_1,'b-o');
set(handles,'LineWidth',1.5);
legend(handles,'分布式干扰对齐');
grid on;
xlabel('SNR(dB)');ylabel('sum rate(bits/s/hz)');
