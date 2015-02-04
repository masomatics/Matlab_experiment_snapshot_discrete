%%
%Preamble
D = genpath('../snap1_26_goodfiles');
addpath(D) 
analysis_load_data
[num_particles, num_all_slices] = size(datx);


%% Construct Compressed snapshots

num_bins = 10; 
field1 = 'weights' ; value1 = zeros(num_bins, num_slices);
field2 = 'values' ; value2 = zeros(num_bins, num_slices);
compress_all_snapshot_hist = struct(field1, value1, field2, value2);
for(k = 1: num_all_slices) 
[num_el, center]= hist(datx(:, k), num_bins);
compress_all_snapshot_hist.weights(:,k) = num_el / sum(num_el);
compress_all_snapshot_hist.values(:,k) = center;
end 
compress_snapshot_hist = ...
    struct(field1, compress_all_snapshot_hist.weights(:,timesample),...
    field2,compress_all_snapshot_hist.values(:,timesample))

%% 
num_iter = 10000;
historylength = 50;
energyhistory = zeros(1,historylength);

theta_now = theta_init;
eta = 0.01;
for(iter = 1 :num_iter)
    rnsource = randn(T,N);
    %Energy E_q[log p] must be getting smaller. 
    [derivative, energy, p_ymk_all] = debug_analysis_mc_deriv(init, ...
theta_now, sigV, sigW, T, rnsource, obs_noise', compress_snap_vals, compress_snap_wgts,...
snaptime, N);
    energyhistory(2:historylength) = energyhistory(1:historylength-1);
    energyhistory(1) = sum(energy);
    theta_now = theta_now + derivative'*eta;
    theta_now(2) = max(min(theta_now(2), 1), -1);
    if(mod(iter, 50)== 0)
        display(['iteration', num2str(iter), ' complete', 'theta =',...
            num2str(theta_now)]);
        display(['energy = ', num2str(sum(energy))]);
        delta_energy = (max(energyhistory) - min(energyhistory))/max(energyhistory)
    end     
    if(delta_energy < 0.1)
    compress_all_snapshot_hist = analysis_make_hists(datx, 50, ...
        num_slices, num_all_slices,timesample)
    break;
    end 
end 

if(delta_energy < 0.1)
    compress_all_snapshot_hist = analysis_make_hists(datx, 50, ...
        num_slices, num_all_slices,timesample)
end 


for(iter = 1 :num_iter)
    rnsource = randn(T,N);
    %Energy E_q[log p] must be getting smaller. 
    [derivative, energy, p_ymk_all] = debug_analysis_mc_deriv(init, ...
theta_now, sigV, sigW, T, rnsource, obs_noise', compress_snap_vals, compress_snap_wgts,...
snaptime, N);
    energyhistory(2:historylength) = energyhistory(1:historylength-1);
    energyhistory(1) = sum(energy);
    theta_now = theta_now + derivative'*eta;
    theta_now(2) = max(min(theta_now(2), 1), -1);
    if(mod(iter, 50)== 0)
        display(['iteration', num2str(iter), ' complete', 'theta =',...
            num2str(theta_now)]);
        display(['energy = ', num2str(sum(energy))]);
        delta_energy = (max(energyhistory) - min(energyhistory))/max(energyhistory)
    end     
end 
  