function [samples_Mt,samples_M] = generate_bernoulli_samples(M,Mt,p,Q)
% 输入 M, Mt, p, Q
% M  每个用户选取天线数
% Mt 每个用户可供选择天线数
% p  概率矩阵
% Q  采样数
% 输出 样本 sample

len = size(p,2);
K = len/Mt; % 用户数
samples_Mt = zeros(Q,len);
samples_M = zeros(Q,M*K);
trial = ones(1,len);
for i = 1:Q
    for j = 1:K
        tmp = binornd(trial(1+(j-1)*Mt:j*Mt),p(1+(j-1)*Mt:j*Mt));
        if sum(tmp) == M     % 恰好满足条件
            samples_Mt(i,1+(j-1)*Mt:j*Mt) = tmp;
        elseif sum(tmp) > M  % 需要随机删除1
            ind = sort(randperm(sum(tmp),sum(tmp)-M));
            num = 0;
            for k = 1:Mt
                if tmp(k) == 1 
                    num = num + 1;
                    if(find(ind == num))
                        tmp(k) = 0;  
                    end
                end
            end
            samples_Mt(i,1+(j-1)*Mt:j*Mt) = tmp;
        else                 % 需要随机添加1
           ind = sort(randperm(Mt-sum(tmp),M-sum(tmp)));
           num = 0;
           for k = 1:Mt
                if tmp(k) == 0 
                    num = num + 1;
                    if(find(ind == num))
                        tmp(k) = 1;  
                    end
                end
           end
           samples_Mt(i,1+(j-1)*Mt:j*Mt) = tmp;
        end
    end
    for j = 1:K
        samples_M(i,1+(j-1)*M:j*M) = find(samples_Mt(i,1+(j-1)*Mt:j*Mt)== 1);
    end
end