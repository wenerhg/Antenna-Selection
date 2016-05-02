% rate_K_user_MIMO(U, H, V, xigma) computes the data rate of a real K-user IC. 
% The inputs are the K zeroforcers U, the K^2 MrxMt channels H, and the K
% beamformers V.

function R = rate_K_user_MIMO(U, H, V,xigma)

R = 0;
[~, d,K] = size(V); % obtain the parameters
S = zeros(d,d,K); % initialize the useful signal matrix
J = zeros(d, (K-1)*d, K); % initialize the interference matrix
  for i = 1:K
        %%% construct interference matrix %%%
        for j = [1:i-1 i+1:K]
             J(:,(j-1)*d+1:j*d,i) = U(:,:,j)'*H(:,:,i,j)*V(:,:,i); %each term is one interference block
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S(:,:,i) = U(:,:,i)'*H(:,:,i,i)*V(:,:,i); % construct the signal space matrix    
        R = R + log2(det(eye(size(J,1))+(1/xigma*U(:,:,i)'*U(:,:,i)+J(:,:,i)*J(:,:,i)')^-1*S(:,:,i)*S(:,:,i)')); % calculate the rate
  end
end

