%%  The  REAL implemetation test.

% This is a beta version of the program that attempts to infer the 
% parameters for the model when the snapshot distribution (empirical)
% q is given. The algorithm follows 5.5 of my snapnote. 
% The theory is given in 2.2. 
%
%
% The model of interest: 
%
% x(m) = theta(1) + theta(2)*x(m-1) + ...
%         theta(3)*x(m-1)/ (1 + x(m-1)^2) + ...
%         theta(4)*cos(theta(5)*x(m-1)) + N(0, sigV^2)
%
% This set of debuggers run the algorithm and creates the movies of 
%     histograms
%     evolving ptilde on snapshots. The heat should be almostall  blue(not much
%     divergence from all 1/num_particles empirical distribution) when the
%     theta00 = theta_a. The heat at point y_0 describes the value 
%     tilde_py(y_0).
% 
%     In particular I note that 
%
%     Nsnap : Number of particles per each snapshot
%     Ntry : Number of  MonteCarlo Xs to be simulated 
%     num_iter : Number of Iterations for Gradient Ascent.

%Optional. Re_extracts the snapshots.
theta_a = [0, 1/2, 25, 6, 0.2];
theta_b = 20;
num_timepts = 40;
sigV = sqrt(10);
%sigW = sqrt(1);
sigW = sqrt(0.4);
Nsnap = 4000; %Number of particles par each snapshot. 
init = -1.3

rnsource_snap = randn(Nsnap, num_timepts);
rnsourcey_snap = randn(Nsnap, num_timepts);
  %Generate data
[xdatm, ydatm]= data_generation_N(init, theta_a, num_timepts, sigV, sigW,  ...
    rnsource_snap, rnsourcey_snap , Nsnap );
timesample = [1, 10:5:40];
  %Generate Snapshots
snapshots = ydatm(:,timesample);


%%
%save('snapshots_real_minus13_0_0_0_0_uneven_aug.mat','snapshots');
hoge = load('snapshots_real_minus13_0_0_0_0_uneven_aug.mat') 
snapshots = hoge.snapshots;

%Initialization. theta_a is the de-facto parameters.
%Set the number of iterations here.
close all;
%theta00 = [0, 1/2, 25, 6, 0.2];
%theta00 = [-0.1327    0.9516    1.09    0.3865    0.2000] 
%theta00 =[ -0.14443      0.5     1.0973     0.36525         0.2]
%theta00 - [-0.0546    0.8765    6.1200    3.1286    1.0426];
%theta00 = [-1, 1/2, 50, 1,  1];

%theta00 = [-1.3,0,0,0,1]; % part 1
%theta00 = [-0.0256    0.8920    8.9758    1.9773    1.0170]; %part2
%theta00 = [-0.1150    0.8308    7.8343    3.9540    1.3202]; %part3
%theta00 = [-0.1301    0.8111    8.4956    4.4360    1.2845]; %part4
%theta00 = [-0.1292    0.7730   10.8109    4.8983    3.9796] %part5

%theta00 = [-0.1934    0.7563   12.8268    5.3425    2.2087]; %part6
%theta00 = [-0.2105    0.7414   13.7000    5.6004    1.3613] %part7
%theta00 = [-0.193174     0.715688      14.1746      5.51424      2.51438] %part8
%theta00=[  -0.1525    0.6351   17.3361    6.1000    7.2770]; %part9
%theta00 = [ -0.1298    0.6773   17.4532    6.0449    8.0412] ; part 10
%theta00 = [   -0.2220    0.6464   18.5362    6.0521    1.4033]; part 11 
%theta00 = [-0.1032    0.6515   18.5006    5.9595   15.5792];
theta00  = [-0.1565    0.6595   19.2218    6.1115    5.1477];
%init = -1.3

num_iter = 3200; 
eta = [0.05, 0.1, 2, 1, 0.2];
[num_particles, num_frames] = size(snapshots);
Ntry = 4000; %Number of montecarlo simulated Xs. 


%Tracks the trajectory of the parameters. 
thetahistory = zeros(num_iter,5);
thetahistory(1,:) = theta00;

%Video Writer  Taking pictures of tilde p and histograms.
writerObj = VideoWriter('snapshot_analysis_Doucet_minus13_0_0_0_0_uneven_part13_aug.avi'); 
writerObj_histo = VideoWriter('snapshot_analysis_Doucet_histogram_minus13_0_0_0_0_uneven_part13_aug.avi'); 

writerObj.FrameRate =2;
writerObj_histo.FrameRate =1;

open(writerObj)
open(writerObj_histo)

snaphistory = zeros(Ntry, num_frames, num_iter);
energyhistory = zeros(1, num_iter) ;

 
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
    energy = zeros(1,length(timesample)-1);
    for(frame_index = 2:length(timesample))
        subplot(2, ceil((length(timesample)-1)/2),frame_index-1)
        hoge = ones(num_particles,1) ;
        %scatter(snapshots(:,frame_index), hoge, 50,  tilde_pys(:,frame_index), 's');
        scatter(snapshots(:,frame_index), tilde_pys(:,frame_index), 50)
        
        energy(frame_index) = mean(log(tilde_pys(:,frame_index)));
        title(['t=', num2str(timesample(frame_index)),'energy = ', num2str(energy(frame_index) ) ]);         
    end
    suptitle(['Heat of tilde_py : Itrn ', num2str(iter),...
       ', Theta=[',num2str(theta00), ']      Total energy =',  num2str(sum(energy))]);
    energyhistory(iter) = sum(energy); 
    
    
    toc
    M(iter) = getframe(gcf);
    writeVideo(writerObj,M(iter));          
    

    %Update theta
    theta00 = theta00 + eta.* deriv_a';
    
    %%%% This sytem is difficult...
            %theta00(3) = ;
            %theta00(4) = 6;
            %theta00(5) = 0.2;
            theta00(2) = min(theta00(2), 1);
            theta00(2) = max(theta00(2), -1);

    %Record.
    thetahistory(iter+1,:) = theta00;
    snaphistory(:,:, iter) = datmat(:, timesample) ;
    
    
    
    figure(2);
    set(gcf,'Position',get(0,'ScreenSize'))

    for frame = 2:length(timesample)
      subplot(2,ceil((length(timesample)-1)/2), frame-1)
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
    display(['The correct parameter is =' , num2str(theta_a)]);

    
end
close(writerObj);
close(writerObj_histo);


save('theta_history_minus13_0_0_0_0_uneven_part13_aug.mat', 'thetahistory')
save('snap_history_minus13_0_0_0_0_uneven_part13_aug.mat', 'snaphistory')
save('energy_history_minus13_0_0_0_0_uneven_part13_aug.mat', 'energyhistory')
