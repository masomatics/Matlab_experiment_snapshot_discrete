function [] = analysis_compare(theta_truth, theta_false, init, num_timepts, sigV, sigW, sys_noise, obs_noise, timesample, N)


[datx_truth, daty_truth]= analysis_data_generation_N(init, theta_truth, num_timepts, sigV, sigW,  sys_noise', obs_noise' , N );
[datx_false, daty_false]= analysis_data_generation_N(init, theta_false, num_timepts, sigV, sigW,  sys_noise', obs_noise' , N );

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

