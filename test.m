clc
clear all
close all
num = 0;
for i = 1:2000
    U = randn(3,3,3);
    V = randn(3,3,3);
    H = randn(3,3,3,3);
    if min_stream_SINR(U,H,V,1) < 0
        num = num + 1;
    end
end
num