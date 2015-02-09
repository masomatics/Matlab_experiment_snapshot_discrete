%% Prepare the histogram data
D = genpath('../snap1_26_goodfiles');
addpath(D) 
analysis_load_data

%% Initialize the paramters 
theta_init = [-1.3, 0,0,0,0.1]; 
%theta = [ 0    0.5000   25.0000    6.0000 0.2000 ]
compress_snapshot_hist = analysis_hist_generation(snapshots, 10);
compress_snap_vals = compress_snapshot_hist.values;
compress_snap_wgts = compress_snapshot_hist.weights;
snaptime = timesample;
T= 40;
N= 15000;

rnsource = randn(T,N);

%%
num_iter = 15000;
theta_now = theta_init;
eta = 0.01;
tic,
for(iter = 1 :num_iter)
    rnsource = randn(T,N);
    %Energy E_q[log p] must be getting smaller. 
    [derivative, energy, p_ymk_all] = analysis_mc_deriv_Doucet_totVar(init, ...
theta_now, sigV, sigW, T, rnsource, compress_snap_vals, compress_snap_wgts,...
snaptime, N);
    theta_now = theta_now - derivative'*eta;
    %theta_now(5) = 0.2;
    theta_now(2) = max(min(theta_now(2), 1), -1);
    if(mod(iter, 100)== 0)
        display(['iteration', num2str(iter), ' complete', 'theta =',...
            num2str(theta_now)]);
        toc,
        display(['energy = ', num2str(sum(energy))]);
        tic,
    end     
end 

toc
%%
analysis_compare(theta, theta_now, init,...
    num_timepts, sigV, sigW, rnsource', rnsource', T, timesample, N, 10)