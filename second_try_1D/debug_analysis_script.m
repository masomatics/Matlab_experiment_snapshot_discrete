theta_init = [-1.3, 0,0,0,0.1]; 
%theta = [ 0    0.5000   25.0000    6.0000 0.2000 ]
compress_snap_vals = compress_snapshot_hist.values;
compress_snap_wgts = compress_snapshot_hist.weights;
snaptime = timesample;
T= 40;
N= 4000;


%%Compute the true energy




%% This will give error when sys_noise, obs_noise is not defined in other piece of debug session
tic, 
theta_now = theta;
[derivative, energy] = debug_analysis_mc_deriv(init, ...
theta, sigV, sigW, T,  sys_noise', obs_noise', compress_snap_vals, compress_snap_wgts,...
snaptime, N)
sum(energy)
toc
%% This will give error when rnsource is not defined in other piece of debug session
tic, 
theta_now = theta;
[derivative, energy] = debug_analysis_mc_deriv(init, ...
theta, sigV, sigW, T, rnsource, rnsource, compress_snap_vals, compress_snap_wgts,...
snaptime, N)
sum(energy)
toc
%%
num_iter = 1000;
theta_now = theta_init;
eta = 0.01;
for(iter = 1 :num_iter)
    rnsource = randn(T,N);
    %Energy E_q[log p] must be getting smaller. 
    [derivative, energy, p_ymk_all] = debug_analysis_mc_deriv(init, ...
theta_now, sigV, sigW, T, rnsource, obs_noise', compress_snap_vals, compress_snap_wgts,...
snaptime, N);
    theta_now = theta_now + derivative'*eta;
    theta_now(2) = max(min(theta_now(2), 1), -1);
    if(mod(iter, 10)== 0)
        display(['iteration', num2str(iter), ' complete', 'theta =',...
            num2str(theta_now)]);
        display(['energy = ', num2str(sum(energy))]);
    end     
end 

%%
close all;
figure(99);
for(frame = 1: num_slices)
    subplot(2, ceil(num_slices/2), frame); 
    plot(compress_snapshot_hist.values(:,frame), [compress_snapshot_hist.weights(:,frame),p_ymk_all(:,frame)] ) 
    legend('truth', 'approx')
end


%%
theta_truth = theta;
%theta_false = [-3, -1, -1, -2, 0.06]; 
%theta_false = theta;
theta_false = [-1.4933          -1    -0.16449   -0.068302     0.13437];
theta_false = [-0.89845           1     0.09717     0.61149    0.053411];
%theta_false = [-3, 0.1,0,0,0];
%theta_false = theta_init
analysis_compare(theta_truth, theta_false, init, num_timepts, sigV, sigW, rnsource, rnsource, timesample, N)
