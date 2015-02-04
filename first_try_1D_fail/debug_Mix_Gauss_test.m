%%
p = 1;
M = 3;
N = 10000;
myu_true = rand(p,M);%[3 1 1];%
sigma_true = rand(p,M);%[2 1 1];%
weight_true = rand(M,1);%[1 1 0];%
weight_true = weight_true/sum(weight_true);
cum_weight = cumsum(weight_true);

for t = 1:N
    % determine the component index by finding
    % the minimum cum_weight that satisfies cum_weight > rand(1)
    [maxc, index] = max((cum_weight > rand(1))./cum_weight);
    
    % generate gaussian sample
    x(:,t) = diag(sqrt(sigma_true(:,index)))*randn(1)+myu_true(:,index); 
end

y = debug_Mix_Gauss(x, myu_true, sigma_true, weight_true);

myu = rand(p,M);
sigma = rand(p,M);
weight = rand(M,1);
weight = weight/sum(weight);
y_rand = debug_Mix_Gauss(x, myu, sigma, weight);

(y > y_rand)
%%