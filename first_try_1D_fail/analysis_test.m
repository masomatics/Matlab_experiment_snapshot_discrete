close all;
theta1 = [  -1.5000         0         0         0    0.2000];
theta1 = rand(1,5)
[deriv_false, energy_false] = analysis_sensitivity(compress_snapshots, timesample,...
    init,theta1, 40, sigV,  sys_noise, 6000)

sum(energy_false)

close all;
theta_true = [0    0.5000   25.0000    6.0000    0.2000];
%theta_true = [-1.5, 0.1, 0,0,0] 
[deriv_true, energy_true] = analysis_sensitivity(compress_snapshots, timesample,...
    init,theta_true, 40, sigV,  sys_noise, 6000)

sum(energy_true)




