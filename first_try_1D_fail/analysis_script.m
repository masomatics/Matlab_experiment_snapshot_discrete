
T_exp = 40;
N_exp = 6000;
learn_rate = 0.1;
num_iter = 10;
theta_init = [-1.5, 0,0,0,0.2];
theta_init = [0, 1/2, 25, 6, 0.2];

theta_now = theta_init;
record = zeros(length(theta_init),num_iter);

tic,

block = 1; 
history_now = 0;
theta_history = zeros(ceil(num_iter/ block), 5); 
energy_history = zeros(ceil(num_iter/ block),1); 



for(iter = 1 : num_iter)
 
    sys_noise2 = randn(N_exp,T_exp);
    [deriv_vals,energy] = analysis_sensitivity(compress_snapshots, timesample,...
        init, theta_now, T_exp, sigV,  sys_noise2, N_exp);
    
    theta_now = theta_now + learn_rate* mean(deriv_vals(2:end,:),1);
    theta_now(2) = min(theta_now(2), 1);
    theta_now(2) = max(theta_now(2), -1);
    theta_now(5) = 0.2; 
    
    if(mod(iter, block) == 0)        
        history_now = history_now + 1
        display(theta_now)
        display(energy)
        theta_history(history_now, :) = theta_now;
        energy_history(history_now) = sum(energy);

    toc
    tic,
    end
    
end
toc
