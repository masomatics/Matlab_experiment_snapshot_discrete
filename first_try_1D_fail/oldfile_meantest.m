
%% Extracts the data.
clear all;
theta_a = [0, 1/2, 25, 6, 0.2];
theta_b = 20;


theta_ax =   [-0.0199    0.8623    7.0718   -3.5525    1.0047]; 
%theta_axx = [-1.3748    -0.84256     -7.6122       2.898      1.3902]
theta_ax = [-0.1150    0.8308    7.8343    3.9540    1.3202]; 
%theta_ax = [-0.0000    0.5108    25.2343    5.9540    0.2002];  % A little perturbation can kill 
%theta_axx = [-0.1301    0.8111    9.4956    4.4360    1.2845];
%theta_axx = [-0.1292    0.7730   10.8109    4.8983    3.9796]
%theta_axx = [-0.1934    0.7563   12.8268    5.3425    2.2087];
%theta_axx =[-0.2105    0.7414   13.7000    5.6004    1.3613];
%theta_axx = [   -0.1932    0.7157   14.1746    5.5142    2.5144]
%theta_axx=[  -0.1525    0.6351   17.3361    6.1000    7.2770]
theta_axx=  [   -0.2220    0.6464   18.5362    6.0521    1.4033];

timesample = [1, 10:5:40];
hoge = load('snapshots_real_minus13_0_0_0_0_uneven.mat')
realmean = mean(hoge.snapshots, 1); 
init = -1.3
%init = -1.3

%theta0 = [-1, 1,10,15,1];
%theta0_b = 20;

sigV = sqrt(10);
sigW = sqrt(1); 



close all;
num_timepts = 40;
Ntry = 100000;


rnsource = randn(Ntry, num_timepts);
rnsource_y = randn(Ntry, num_timepts);
[xdatm, ydatm]= data_generation_N(init, theta_a, num_timepts, sigV, sigW,  rnsource, rnsource_y , Ntry );


[xdatmx, ydatmx]= data_generation_N(init, theta_ax, num_timepts, sigV, sigW,  rnsource, rnsource_y , Ntry );
[xdatmxx, ydatmxx]= data_generation_N(init, theta_axx, num_timepts, sigV, sigW,  rnsource, rnsource_y , Ntry );



for k = 1: num_timepts

ydatmean(k) = mean(ydatm(:,k)); 
ydatmeanx(k) = mean(ydatmx(:,k)); 
ydatmeanxx(k) = mean(ydatmxx(:,k)); 

end

%%
%close all;
timesample = [1, 15:5:40];

figure(99)
    set(gcf,'Position',get(0,'ScreenSize'))
set(gca,'FontSize',20)
plot(1:num_timepts, ydatmean, 'r');
hold on;
plot(1:num_timepts, ydatmeanx, 'g');
hold on;
plot(1:num_timepts, ydatmeanxx, 'b');
hold on;
plot(timesample, realmean, 'bo');

plot(timesample, ydatmean(timesample), 'ko',  'MarkerFaceColor',[.1 0 .0], 'MarkerSize', 7);
yL = get(gca,'YLim');

for m = 1:length(timesample)

line([timesample(m), timesample(m)],yL,'Color','k');

end 
xlabel('n', 'FontSize',20) 
legend('E[Real Dynamics]','E[Esimated Dynamics]', 'Time pts of Snapshots')

hold off;
