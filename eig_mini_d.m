function [varargout]=eig_mini_d(A,d)
%%用于产生d个最小特征值对应的特征向量
[V, D] = eig(A);
D = abs(diag(D));
[D, I] = sort(D, 'ascend');
if d > length(D)
    d = length(D);
end
%varargout = {D(1 : d)};
varargout = {V(:, I(1 : d))};
end