% maxSINR_V(H, V, P_stream, iter) runs the max-SINR algorithm and returns 
% the set of beamforming and zeroforcing matrices V and U. 
% The inputs are the K^2 MrxMt channels H, the K initialization 
% beamformers V, the power allocated per user's stream P_stream, and the number of iterations.

function [V, U] = maxSINR_K_user(H, V, xigma, iter)

[Mr, Mt, K, ~] = size(H); % obtain problem parameters
d = size(V,2);
U = zeros(Mr, d, K);

for i = 1:iter
    for k = 1:K
        for l = 1:d
            B_kl = 0; % initialize interference plus noise covariance matrix
            % compute the interference plus noise covariance matrix
            for j = 1:K
                for d_i = 1:d
                    B_kl = B_kl + xigma*H(:,:,k,j)*V(:,d_i,j)*V(:,d_i,j)'*H(:,:,k,j)';
                end
                 B_kl = B_kl - xigma*H(:,:,k,k)*V(:,l,k)*V(:,l,k)'*H(:,:,k,k)';
            end
           
            B_kl = B_kl+eye(Mr);
            %%%
            B_kl_inv = B_kl^-1;
            U(:,l,k) = B_kl_inv*H(:,:,k,k)*V(:,l,k)./ norm(B_kl_inv*H(:,:,k,k)*V(:,l,k)); % compute the l-th column of the k-th zeroforcer
        end
    end
    for k = 1:K
        for l = 1:d
            B_kl = 0; % initialize reciprocal interference plus noise covariance matrix
            % compute the reciprocal interference plus noise covariance matrix
            for j = 1:K
                for d_i = 1:d
                    B_kl = B_kl + xigma*H(:,:,j,k)'*U(:,d_i,j)*U(:,d_i,j)'*H(:,:,j,k);
                end
                B_kl = B_kl - xigma*H(:,:,k,k)'*U(:,l,k)*U(:,l,k)'*H(:,:,k,k);
            end
             
            %%%
            B_kl = B_kl+eye(Mt);
            B_kl_inv = B_kl^-1;
            V(:,l,k) = B_kl_inv*H(:,:,k,k)'*U(:,l,k)./ norm(B_kl_inv*H(:,:,k,k)'*U(:,l,k)); % compute the l-th column the k-th beamformer
        end
    end   
end

end