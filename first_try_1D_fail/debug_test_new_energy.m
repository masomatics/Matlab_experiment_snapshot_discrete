%%
N = 10000
q = randn(1,N).^2;
q = q/sum(q);
p = randn(1,N).^2;
p = p/sum(p);
dim = 1

pp = [q == max(q)];
sum(q.* log(q) )- log(sum(q.*q)) > sum(p.* log(q) )- log(sum(p.*q))
sum(q.* log(q) )- log(sum(q.*q)) > sum(pp.* log(q) )- log(sum(pp.*q))


%%
analysis_load_data
%%
[num_par,M, num_slice] = size(compress_snapshots)
myu = rand(dim,M);
sigma = rand(dim,M);
weight = rand(M,1);
weight = weight/sum(weight);
cum_weight = cumsum(weight_true);


for t = 1:N
    % determine the component index by finding
    % the minimum cum_weight that satisfies cum_weight > rand(1)
    [maxc, index] = max((cum_weight > rand(1))./cum_weight);
    
    % generate gaussian sample
    x(:,t) = diag(sqrt(sigma_true(:,index)))*randn(1)+myu_true(:,index); 
end


y = debug_Mix_Gauss(x, myu_true, sigma_true, weight_true)
objectiv_truth = 0; 
 for(m = 1 : num_components)
            objectiv_truth = objectiv_truth + GMMtruth(m, 3) * ...
                exp(-(x - GMMtruth(m, 1) ).^2/ ...
                (2*GMMtruth(m, 2))...
                )./sqrt(GMMtruth(m, 2)); 
 end

 
 
 
 
%xfake = 1.5*randn(1,N);
xfake = debug_GMM_generate(GMMfake1,N); 
yfake = debug_Mix_Gauss(xfake, myu_true, sigma_true, weight_true);
objectiv_fake = 0; 
 for(m = 1 : num_components)
            objectiv_fake = objectiv_fake + GMMtruth(m, 3) * ...
                exp(-(xfake - GMMtruth(m, 1) ).^2/ ...
                (2*GMMtruth(m, 2))...
                )./sqrt(GMMtruth(m, 2)); 
 end
