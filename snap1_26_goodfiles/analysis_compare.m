% Compares two Models with different theta for the model 
%
% The model of interest: 
%
% x(m) = theta(1) + theta(2)*x(m-1) + ...
%         theta(3)*x(m-1)/ (1 + x(m-1)^2) + ...
%         theta(4)*cos(theta(5)*x(m-1)) + N(0, sigV^2)
% y(m) = N(x(m), sigW^2);
% x(1) = N(initmean, sigW^2)
% Inputs: 
%   initmean    :As described above
%   theta_true       :First Choice of the theta. Used as q in E[log p]
%   theta_false      :Second Choice of the theta
%   sigV        :System noise.
%   sigW        :Observation noise. 
%   sys_noise   : Standard Normal (N, timesample) 
%   obs_noise   : Standard Normal (N, timesample)
%   T           :Terminal time. Integer
%   timesample  :times at which the snapshots were taken. The length is
%               num_frames.
%   N           :number of Xpaths to be sampled
%   nbins       :Value used for histogram comparison
%
% Outputs: 
%               JUST PLOTS
%
%               Plot97 : On the support of the hist(true snapshot,
%               nbins),create the distribution on p from false snapshot.
%               tilde p in the SnapNote.
%               Plot98 : Two histograms covered on the same graph
%               Plot99 : Moment comparison

function [] = analysis_compare(theta_truth, theta_false, init,...
    num_timepts, sigV, sigW, sys_noise, obs_noise, T, timesample, N, nbins)

close all;
[datx_truth, daty_truth]= analysis_data_generation_N(init, theta_truth, num_timepts, sigV, sigW,  sys_noise, obs_noise , N );
[datx_false, daty_false]= analysis_data_generation_N(init, theta_false, num_timepts, sigV, sigW,  sys_noise, obs_noise , N );

hoge = figure(98);
num_slices = length(timesample)

snapshots_truth = daty_truth(:,timesample);
snapshots_false = daty_false(:,timesample);

figure(98)
    for(frame = 1: num_slices)
        subplot(2, ceil(num_slices/2), frame)
        num_breaks = 50;
        h1 = histogram(snapshots_truth(:,frame),num_breaks);
        hold on;
        h2 = histogram(snapshots_false(:,frame),num_breaks);    
        title(['t=', num2str(timesample(frame))])
        legend('truth','false')
    end

m1_data_true =  mean(snapshots_truth, 1);
m2_data_true = mean(snapshots_truth.^2, 1);
var_data_true= var(snapshots_truth, 1);
m3_data_true = mean(snapshots_truth.^3, 1);
m4_data_true = mean(snapshots_truth.^4, 1);
m5_data_true = mean(snapshots_truth.^5, 1);
m6_data_true = mean(snapshots_truth.^6, 1);
    
m1_data_approx =  mean(snapshots_false, 1);
m2_data_approx = mean(snapshots_false.^2, 1);
var_data_approx= var(snapshots_false, 1);
m3_data_approx = mean(snapshots_false.^3, 1);
m4_data_approx = mean(snapshots_false.^4, 1);
m5_data_approx = mean(snapshots_false.^5, 1);
m6_data_approx = mean(snapshots_false.^6, 1);

figure(99)
subplot(2,3,1)
plot(timesample, m1_data_true)
hold on;
plot(timesample, m1_data_approx)
title('mean')


subplot(2,3,2)

plot(timesample, var_data_true)
hold on;
plot(timesample, var_data_approx)
title('var')


subplot(2,3,3)

plot(timesample, m3_data_true)
hold on;
plot(timesample, m3_data_approx)
title('moment3')



subplot(2,3,4)

plot(timesample, m4_data_true)
hold on;
plot(timesample, m4_data_approx)
title('moment4')


subplot(2,3,5)

plot(timesample, m5_data_true)
hold on;
plot(timesample, m5_data_approx)
title('moment5')


subplot(2,3,6)

plot(timesample, m6_data_true)
hold on;
plot(timesample, m6_data_approx)
title('moment6')

figure(97)
compress_snapshot_hist = analysis_hist_generation(snapshots_truth, nbins)
compress_snap_vals= compress_snapshot_hist.values;
compress_snap_wgts= compress_snapshot_hist.weights;
[derivative, energy, p_ymk_all] = analysis_mc_deriv_Doucet(init, ...
theta_false, sigV, sigW, T, sys_noise', compress_snap_vals, compress_snap_wgts, timesample, N);

    for(frame = 1: num_slices)
        subplot(2, ceil(num_slices/2), frame)
        plot(compress_snap_vals(:,frame), [compress_snap_wgts(:,frame), p_ymk_all(:,frame)]) 
        legend('truth','false') 
        energy_true = log(compress_snap_wgts(:,frame))'*compress_snap_wgts(:,frame);
        energy_false= energy(frame);
        title(['TrueE', num2str(energy_true), ' FalseE=', num2str(energy_false)]) 
    end

end 

