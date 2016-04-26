% normalize(X) outputs the normalized matrices along the 3rd dimension of a 3-D
% matrix.
% The input is the 3-D matrix X.

function Xnorm = normalize(X)

Xnorm = zeros(size(X));
for i = 1:size(X,3)
    Xnorm(:,:,i) = X(:,:,i)./norm(X(:,:,i)); % normalize X(:,:,i)
end

end

