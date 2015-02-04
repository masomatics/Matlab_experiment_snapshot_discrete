%% Prepare the histogram data
analysis_hist_snapshots
%% Initialize the paramters 
theta_init = [-1.3, 0,0,0,0.1]; 
%theta = [ 0    0.5000   25.0000    6.0000 0.2000 ]
compress_snap_vals = compress_snapshot_hist.values;
compress_snap_wgts = compress_snapshot_hist.weights;
snaptime = timesample;
T= 40;
N= 3000;

rnsource = randn(T,N);

%% The derivative at the Real parameter
tic, 
derivative = analysis_mc_deriv(init, ...
theta, sigV, sigW, T, rnsource, compress_snap_vals, compress_snap_wgts,...
snaptime, N)
toc
tic,
derivative = oldfile_analysis_mc_deriv(init, ...
theta, sigV, sigW, T, rnsource, compress_snap_vals, compress_snap_wgts,...
snaptime, N)
toc
%% Iteration of seepest ascent
num_iter = 1000;
theta_now = theta_init;
eta = 0.01;
for(iter = 1 :num_iter)
    derivative = analysis_mc_deriv(init, ...
theta_now, sigV, sigW, T, rnsource, compress_snap_vals, compress_snap_wgts,...
snaptime, N);
    theta_now = theta_now + derivative'*eta;
    theta_now(2) = max(min(theta_now(2), 1), -1);
    if(mod(iter, 10)== 0)
        display(['iteration', num2str(iter), ' complete', 'theta =',...
            num2str(theta_now)]);
    end     
end 