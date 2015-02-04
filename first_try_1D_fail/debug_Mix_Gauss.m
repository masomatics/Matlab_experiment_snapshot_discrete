function y = Mix_Gauss(x, myu, sigma, weight)

% x: p*N matrix, N samples of p-dim observation data
%
% M-components Gaussian mixtures
% myu: p*M matrix, mean of M-components Gaussian mixtures
% sigma: p*p*M (or p*M) matrix, covarianve of M-components Gaussian mixtures
% weight: M*1 vector. weights of M-components Gaussian mixtures

N = size(x,2);
M = length(weight);
f = zeros(N,M);

for m=1:M
    switch ndims(sigma)
        case 1
            inv_sigma = 1/sigma(m);
        case 2
            inv_sigma = diag(1./squeeze(sigma(:,m)));
        case 3
            inv_sigma = inv(squeeze(sigma(:,:,m)));
    end
    
    temp = (x-repmat(myu(:,m),1,N));
    temp2 = temp'*inv_sigma;
    for n = 1:N
        inner_prod(n) = temp2(n,:)*temp(:,n);
    end
    f(:,m) = exp(-inner_prod(:)/2)*sqrt(det(inv_sigma));
end

y = mean(log(f*weight(:)));