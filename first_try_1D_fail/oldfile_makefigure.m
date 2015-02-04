clear all
timepicks = [1,100, 200] ;
hoge = load('thetahistory_initial_1_1_1_1_02.mat');
thetahistory_p1 = hoge.thetahistory;

%Snapshot extraction
theta_a = [0, 1/2, 25, 6, 0.2];
theta_b = 20;
num_timepts = 40;
sigV = sqrt(10);
sigW = sqrt(0.4);
Nsnap = 4000; %Number of particles par each snapshot. 
Ntry = 4000;
init = -1.5

rnsource_snap = randn(Nsnap, num_timepts);
rnsourcey_snap = randn(Nsnap, num_timepts);
  %Generate data
[xdatm, ydatm]= data_generation_N(init, theta_a, num_timepts, sigV, sigW,  ...
    rnsource_snap, rnsourcey_snap , Nsnap );
timesample = [1, 5:5:40];
  %Generate Snapshots
snapshots = ydatm(:,timesample);

close all;
for m = 1 : length(timepicks)
    
    theta00 = thetahistory_p1(timepicks(m),:) 
    Xf = figure(1+2*(m-1));
    set(gcf,'Position',get(0,'ScreenSize'))

    %Genenrate the normal rancdom variables
    rnsource = randn(Ntry, num_timepts);
    
    %The derivative code
    [datmat,tilde_pys, deriv_a] = data_generation_N_deriv_stat_beta...
      (init, theta00,theta_b, sigV, sigW, num_timepts, rnsource, snapshots, timesample, Ntry);  
   
    %Plots the tilde_p on the snapshots by heat diagram.
    energy = zeros(1,length(timesample)-1);
    for(frame_index = 2:length(timesample))
        subplot(2, (length(timesample)-1)/2,frame_index-1)
        %scatter(snapshots(:,frame_index), hoge, 50,  tilde_pys(:,frame_index), 's');
        scatter(snapshots(:,frame_index), tilde_pys(:,frame_index), 50)
        
        energy(frame_index) = mean(log(tilde_pys(:,frame_index)));
        title(['t=', num2str(timesample(frame_index)),'energy = ', num2str(energy(frame_index) ) ]);         
    end
    suptitle(['Heat of tilde_py : Itrn ', num2str(timepicks(m)),...
       ', Theta=[',num2str(theta00), ']      Total energy =',  num2str(sum(energy))]);
   
       snaprecord(:,:, m) = datmat(:, timesample) ;
       filename = ['Itr', num2str(timepicks(m)),'histogram111102.jpg']; 
    saveas(Xf,filename);
       
       
    Xf2 = figure(2*m);
    set(gcf,'Position',get(0,'ScreenSize'))

    for frame = 2:9
      subplot(2,4, frame-1)
      hist(snapshots(:,frame),100);
      h = findobj(gca,'Type','patch');
      set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
      hold on;
      hist(snaprecord(:,frame,m), 100); 
      h = findobj(gca,'Type','patch');
      set(h,'facealpha',0.75);
      %set(h(1),'FaceColor','r','EdgeColor','k');
      %set(h(2),'FaceColor','b','EdgeColor','k');
      % set(h(2),'FaceColor','b','EdgeColor','k');
      title(['t=', num2str(timesample(frame))]);
      legend('Observed','Monte Carlo')
      hold off;
    end
        suptitle(['Heat of tilde_py : Itrn ', num2str(timepicks(m)),...
       ', Theta=[',num2str(theta00), ']      Total energy =',  num2str(sum(energy))]);
   
          filename2 = ['Itr', num2str(timepicks(m)),'histogram111102_energy.jpg']; 
    saveas(Xf2,filename2);
    
end 