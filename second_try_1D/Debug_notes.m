%%
NN = 3;
rnsource3 = rand(T,NN) ;
[derivative, energy] = debug_analysis_mc_deriv(init, ...
theta, sigV, sigW, 40, rnsource3, rnsource3, compress_snap_vals, compress_snap_wgts,...
snaptime, NN)


%%
rnsource = sys_noise';
rnsource2= obs_noise';
compress_snap_vals = compress_snapshot_hist.values;
compress_snap_wgts = compress_snapshot_hist.weights;
snaptime = timesample;

[derivative, energy] = debug_analysis_mc_deriv(init, ...
theta, sigV, sigW, 40, rnsource, rnsource2, compress_snap_vals, compress_snap_wgts,...
snaptime, N)

plot(compress_snap_vals(:,snaptime_now),compress_snap_wgts(:,snaptime_now)*N)
hold on;
%[hoge, piyo] = hist(daty);
[hoge, piyo] =hist(datmat);
plot(piyo, hoge)
[hoge, piyo] =hist(datmat + w(time+1,:));
plot(piyo, hoge)
[hoge, piyo] =hist(datmat + w(time-1,:));
plot(piyo, hoge)
hold off;


%% Validating the derivative
N = 5000000
rnsource_debug = randn(T,N);
theta_now = theta_init
%theta_now = theta_init
epsilon = 0.01
tic,
[derivative, energy] = debug_analysis_mc_deriv(init, ...
theta_now, sigV, sigW, T, rnsource_debug, rnsource_debug, compress_snap_vals, compress_snap_wgts,...
snaptime, N)
toc
tic,
[derivative, energy] = debug_analysis_mc_deriv_temp(init, ...
theta_now, sigV, sigW, T, rnsource_debug, rnsource_debug, compress_snap_vals, compress_snap_wgts,...
snaptime, N)
toc

tic,
derivative = analysis_mc_deriv(init, ...
theta_now, sigV, sigW, T, rnsource_debug, compress_snap_vals, compress_snap_wgts, snaptime, N)
toc

% tic,
% [derivative, energy] = oldfile_debug_analysis_mc_deriv_minus1(init, ...
% theta_now, sigV, sigW, T, rnsource_debug,rnsource_debug, compress_snap_vals, compress_snap_wgts, snaptime, N)
% toc
%%
for(k = 1 :5) 
    theta_new = theta_now;
    theta_new(k) = theta_new(k) + epsilon;
    tic,
    [derivative_new, energy_new] = debug_analysis_mc_deriv_temp(init, ...
theta_new, sigV, sigW, T, rnsource_debug, rnsource_debug, compress_snap_vals, compress_snap_wgts,...
snaptime, N);

%     [derivative_new, energy_new] = oldfile_debug_analysis_mc_deriv_minus1(init, ...
% theta_new, sigV, sigW, T, rnsource_debug, rnsource_debug, compress_snap_vals, compress_snap_wgts,...
% snaptime, N);
%     toc,
    display(['compare ',num2str([(sum(energy_new) - sum(energy))/epsilon, derivative(k)])]);
end 


%% Validating the derivative II 
N = 5000000
rnsource_debug = randn(T,N);
theta_now = theta
tic,
[derivative,deriv_all, statmean] = debug_analysis_mc_deriv_simplestat(init, ...
theta, sigV, T,  rnsource_debug,timesample, N)
toc

realstatmean =  mean(datx, 1);
statmean(1) = init;
%plot(timesample, [realstatmean(timesample); statmean])
%
epsilon = 0.01;
N = 5000000
rnsource_debug = randn(T,N);
for(k = 1 :5) 
    theta_new = theta_now;
    theta_new(k) = theta_new(k) + epsilon;
    [derivative_new, deriv_all_new statmean_new] = debug_analysis_mc_deriv_simplestat(init, ...
theta_new, sigV, T,  rnsource_debug,timesample, N);
    statmean_new(1) = init;
    derivative_fd(k) = (sum(statmean_new) - sum(statmean))/epsilon;
    derivative_all_fd(k,:) = (statmean_new - statmean)/epsilon;
end 
    display([ derivative_fd; derivative']);

    
    
%%     
    