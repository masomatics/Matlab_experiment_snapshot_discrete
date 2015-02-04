%% Extracts the data.
theta_a = [0, 1/2, 25, 6, 0.2];
theta_b = 20;


init = -1.5

%theta0 = [-1, 1,10,15,1];
%theta0_b = 20;

sigV = sqrt(10);
sigW = sqrt(1); 


datx= load('xdata6000particles.mat');
daty= load('ydata6000particles.mat');

datx= datx.xdatm;
daty= daty.ydatm;

timesample = [1, 5:5:40];

snapshots = daty(:,timesample);



%% This part checks if mean is working correctly. Try with Ntry > 6000.
%The saved data over 6000 can be a little on rough
%side,
close all;
num_timepts = 40;
Ntry = 100000;
%sys_noise = randn(1,num_timepts);
%obs_noise = randn(1,num_timepts);

rnsource = randn(Ntry, num_timepts);
rnsource_y = randn(Ntry, num_timepts);

[xdatm, ydatm]= data_generation_N(init, theta_a, num_timepts, sigV, sigW,  rnsource, rnsource_y , Ntry );

%[datmat, tilde_pys, deriv_a] = data_generation_N_deriv_stat...
%      (init, theta_a,theta_b, sigV, sigW, num_timepts, rnsource, snapshots, timesample, Ntry);

  
%[tilde_pys, deriv_a] = data_generation_N_deriv_stat_beta...
%      (init, theta_a,theta_b, sigV, sigW, num_timepts, rnsource, snapshots, timesample, Ntry);  
  
figure(1);
for k = 1: num_timepts

xdatmean(k) = mean(xdatm(:,k)); 
ydatmean(k) = mean(ydatm(:,k)); 
datm_mean(k) = mean(daty(:,k));

end

plot(timesample, xdatmean(timesample), 'b');
%ylim([min(xdatmean(timesample)-1),max(xdatmean(timesample)+1)])

hold on;
plot(timesample, ydatmean(timesample), 'r');
plot(timesample, datm_mean(timesample),'g');

ylim([min(ydatmean(timesample)-2),1])

snaps = tilde_pys(:,2:end);

hold off; 
figure(5);

writerObj = VideoWriter('snapshots_Y_2.avi');
writerObj.FrameRate = 1;


open(writerObj);
nbins = 300;


for frame =1:40
    hist(ydatm(:,frame), nbins)
        ylim([0,3*N/nbins])
        xlim([ymin,ymax])

    M(frame) = getframe(gcf);
    writeVideo(writerObj,M(frame));
end 
close(writerObj)



%%  The  REAL implemetation test.


%
% Runs the algorithm and creates the movies of 
%     histograms
%     evolving ptilde on snapshots. The heat should be almostall  blue(not much
%     divergence from all 1/num_particles empirical distribution) when the
%     theta00 = theta_a. The heat at point y_0 describes the value 
%     tilde_py(y_0).


%Optional. Re_extracts the snapshots.
theta_a = [0, 1/2, 25, 6, 0.2];
theta_b = 20;
num_timepts = 40;
sigV = sqrt(10);
sigW = sqrt(3);
Nsnap = 3000;
rnsource_snap = randn(Nsnap, num_timepts);
rnsourcey_snap = randn(Nsnap, num_timepts);
  %Generate data
[xdatm, ydatm]= data_generation_N(init, theta_a,num_timepts, sigV, sigW,  ...
    rnsource_snap, rnsourcey_snap , Nsnap );
timesample = [1, 5:5:40];
  %Generate Snapshots
snapshots = ydatm(:,timesample);

%Initialization. theta_a is the de-facto parameters.
%Set the number of iterations here.
close all;
theta00 = [1,1,1,1,0.2];
theta00 = theta_a; 
num_iter = 1; 
eta = 0.1;
[num_particles, num_frames] = size(snapshots);
Ntry = 3000;
%Tracks the trajectory of the parameters. 
thetahistory = zeros(num_iter,5);
thetahistory(1,:) = theta00;

%Video Writer  Taking pictures of tilde p and histograms.
writerObj = VideoWriter('snapshot_analysis_Doucet.avi'); 
writerObj_histo = VideoWriter('snapshot_analysis_Doucet_histogram.avi'); 

writerObj.FrameRate =2;
writerObj_histo.FrameRate =1;

open(writerObj)
open(writerObj_histo)

snaphistory = zeros(Ntry, num_frames, num_iter);
for iter = 1:num_iter
    tic,
    figure(1);
    set(gcf,'Position',get(0,'ScreenSize'))

    %Genenrate the normal rancdom variables
    rnsource = randn(Ntry, num_timepts);
    
    %The derivative code
    [datmat,tilde_pys, deriv_a] = data_generation_N_deriv_stat_beta...
      (init, theta00,theta_b, sigV, sigW, num_timepts, rnsource, snapshots, timesample, Ntry);  
   
    %Plots the tilde_p on the snapshots by heat diagram.
    for(frame_index = 2:length(timesample))
        subplot(1,(length(timesample)-1),frame_index-1)
        hoge = ones(num_particles,1) ;
        scatter(hoge,snapshots(:,frame_index), 8, tilde_pys(:,frame_index));
        title(['t=', num2str(timesample(frame_index))]);
    end
    suptitle(['Heat of tilde_py : Itrn ', num2str(iter), ',   Theta=[',num2str(theta00), ']']);

    
    
    toc
    M(iter) = getframe(gcf);
    writeVideo(writerObj,M(iter));          
    

    %Update theta
    theta00 = theta00 + eta* deriv_a';
    theta00(5) = 0.2;

    %Record.
    thetahistory(iter+1,:) = theta00;
    snaphistory(:,:, iter) = datmat(:, timesample) ;
    
    
    
    figure(2);
    set(gcf,'Position',get(0,'ScreenSize'))

    for frame = 2:9
      subplot(2,4, frame-1)
      hist(snapshots(:,frame),100);
      h = findobj(gca,'Type','patch');
      set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
      hold on;
      hist(snaphistory(:,frame,iter), 100); 
      h = findobj(gca,'Type','patch');
      set(h,'facealpha',0.75);
      %set(h(1),'FaceColor','r','EdgeColor','k');
      %set(h(2),'FaceColor','b','EdgeColor','k');
      % set(h(2),'FaceColor','b','EdgeColor','k');
      title(['t=', num2str(timesample(frame))]);
      legend('Observed','Monte Carlo')
      hold off;
    end
    M_h(iter) = getframe(gcf);
    writeVideo(writerObj_histo,M_h(iter));  
    
    
    
    
    display([num2str(iter), 'th iteration:deriv is =' , num2str(deriv_a')]);

    display([num2str(iter), 'th iteration complete:  new theta0 is =' , num2str(theta00)]);
    
    
end
close(writerObj);
close(writerObj_histo);


%% For debugging purpose. Compares the actual snapshots
%figure(2);
for k = 1: num_timepts
    

datmean(k) = mean(datmat(:,k)); 


end

 % Compare the 
figure(3);
set(gcf,'Position',get(0,'ScreenSize'))
for frame = 2:9
  subplot(2,4, frame-1)
  hist(snapshots(:,frame),100);
  h = findobj(gca,'Type','patch');
  set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
  hold on;
  hist(snaphistory(:,frame,iter), 100); 
  h = findobj(gca,'Type','patch');
  set(h,'facealpha',0.75);
  %set(h(1),'FaceColor','r','EdgeColor','k');
 %set(h(2),'FaceColor','b','EdgeColor','k');
 % set(h(2),'FaceColor','b','EdgeColor','k');
  title(['t=', num2str(timesample(frame))]);

  hold off;
end

figure(4);
plot(1:num_timepts, datmean,'b')
hold on;
plot(timesample, datmean(timesample), 'r')