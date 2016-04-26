function [min_SINR] = min_stream_SINR(U, H, V,xigma)

[Mr, ~, K, ~] = size(H); % obtain problem parameters
d = size(V,2);
stream_SINR = zeros(d,K);



for k = 1:K
    for l = 1:d
        B_kl = 0; % initialize interference plus noise covariance matrix
        % compute the interference plus noise covariance matrix
        for j = 1:K
            for d_i = 1:d
                B_kl = B_kl + xigma*H(:,:,k,j)*V(:,d_i,j)*V(:,d_i,j)'*H(:,:,k,j)';
            end          
        end
        B_kl = B_kl - xigma*H(:,:,k,k)*V(:,l,k)*V(:,l,k)'*H(:,:,k,k)';
        B_kl = B_kl+eye(Mr);
        stream_SINR(l,k) = xigma*U(:,l,k)'*H(:,:,k,k)*V(:,l,k)*V(:,l,k)'*H(:,:,k,k)'*U(:,l,k)./(U(:,l,k)'*B_kl*U(:,l,k));
    end
end
min_SINR = min(stream_SINR(:));