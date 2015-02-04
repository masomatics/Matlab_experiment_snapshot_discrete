% data_generation_N_deriv_stat_beta
%
% This is a beta version of the program that computes the derivative of the
% objective function 
% E_q [log_\theta p(y, theta)]   
% when q is given. The algorithm follows 5.5 of my snapnote. 
% The theory is given in 2.2. 
%
%
% The model of interest: 
%
% x(m) = theta(1) + theta(2)*x(m-1) + ...
%         theta(3)*x(m-1)/ (1 + x(m-1)^2) + ...
%         theta(4)*cos(theta(5)*x(m-1)) + N(0, sigV^2)
%
% x(1) = N(initmean, sigW^2)
%
% Inputs: 
%   initmean    :As described above
%   theta       :As described above
%   theta_b     :Not used in this code.       
%   sigV        :System noise.
%   sigW        :Observation noise. 
%   T           :Terminal time. Discrete.
%   rnsource    :std normal randoms from which the path of x will be generated.
%   snapshots   :snapshots, in matrix. num_particles( of y) by num_frames.
%   snaptime    :times at which the snapshots were taken. The length is
%               num_frames.
%   N           :number of Xpaths to be sampled
%
% Outputs: 
%
%   datmat      :Stores all trajectories of x. N by T. 
%   tilde_pys   :the approx of p in 2.(a).(iv) of pseudocode. num_particles
%                by num_frames.
%   deriv_a     :Nabla theta.  length is 5.

function [datmat, tilde_pys, deriv_a] = data_generation_N_deriv_stat_beta...
 (initmean, theta,theta_b, sigV, sigW, T, rnsource, snapshots, snaptime, N)

  [num_particles, num_frames] = size(snapshots);
  
  default = false;
  if(default == true)
    %N = 100;
    %    a = [1/2, 25, 8, 1.2];
    theta = [0,1/2, 25, 6, 0.2];
    b = 20;
    %    b =10;   
    sigV = sqrt(10);
    sigW = sqrt(3);
    %sigW = sqrt(5);
  end 
  
  %%%Global variables%%%
  %parameters = {theta, b, sigV, sigW};
  parameters = {theta, theta_b};
  display(parameters);
  time =  1:T  ;
  
  
  %%%Target variables%%%% 
  datmat = NaN(N,T);
  deriv_a= zeros(5,1);
  deriv_a_frames= zeros(5,num_frames);  %\nabla at each snapshot.
  tilde_pys = zeros(num_particles, num_frames); %2.(a).iv before the av.
  tilde_pym = zeros(num_particles, N);
  deriv_pym = zeros(num_particles, N, 5); %2.(a).iii before the average
  
  
  %%%Initialization%%%
  v = sigV*rnsource;
  initx = sqrt(5)* randn(N,1)  + initmean;
  %initx = initmean;
  datmat(:,1) = initx;
  snaptime_now = 2;  
  deriv_loglike = zeros(N,5);

  
  %%%Temporary variables%%%
  %These are just holders of intermediary variavles.
  ftheta = NaN(N,5);
  datmathat = initx;
  deriv_logtemp = zeros(num_particles,N,5);
  tilde_pytemp = zeros(num_particles,N,5);
  tilde_pytemp2= zeros(num_particles,1,5);
  temp1 = zeros(num_particles,N);
  temp2 = zeros(num_particles,N);
  temp3 = zeros(num_particles,N);
  
  for(m = 2:T)

    %Updating  the states
    h = waitbar(m/T) ;
    datmathat = theta(1) + theta(2).*datmat(:,m-1) + ...
         theta(3).*datmat(:,m-1)./ (1 + (datmat(:,m-1)).^2) + ...
         theta(4).*cos(theta(5).*datmat(:,m-1)) ;
    
    %Temporary variables to be used for deriv_loglike. The elemsnts of
    % \nabla_theta f(x(m-1))
    ftheta(:,1) = ones(N,1); 
    ftheta(:,2) = datmat(:,m-1);
    ftheta(:,3) = datmat(:,m-1)./ (1 + (datmat(:,m-1)).^2);
    ftheta(:,4) = cos(theta(5).*datmat(:,m-1));
    ftheta(:,5) = theta(4).*datmat(:,m-1).*(-sin(theta(5).*datmat(:,m-1)));
     
    % x(t) = f(x(m-1))+ epsilon 
    datmat(:,m) = datmathat + v(:,m);
    
    % \nabla_theta log p( x[0,t_m] )  
    deriv_loglike = deriv_loglike + ...
        repmat((datmat(:,m) - datmathat)/(sigV^2),1,5).*ftheta;
    
    if(m == snaptime(snaptime_now))
        m ;
%        % tic,
%         for j = 1 : N 
%             for k = 1 : num_particles
%                 yxtemp(k,j) = snapshots(k,snaptime_now) - datmat(j,m) ;
%                 yx(k,j) = min(yxtemp(k,j)^2/(2*sigW^2), 600);
%                 tilde_pym2(k,j)  =  exp(-yx(k,j)); 
%                 deriv_pym2(k,j,1)=  deriv_loglike(j,1)*tilde_pym2(k,j) ; 
%                 deriv_pym2(k,j,2)=  deriv_loglike(j,2)*tilde_pym2(k,j) ;
%                 deriv_pym2(k,j,3)=  deriv_loglike(j,3)*tilde_pym2(k,j) ;
%                 deriv_pym2(k,j,4)=  deriv_loglike(j,4)*tilde_pym2(k,j) ;
%                 deriv_pym2(k,j,5)=  deriv_loglike(j,5)*tilde_pym2(k,j) ;
%             end 
%         sum_tilde_pym2(j)= sum(tilde_pym2(:,j));
%         deriv_pym2(:,j,:) = deriv_pym2(:,j,:)/sum_tilde_pym2(j);
%         tilde_pym2(:,j) = tilde_pym2(:,j)/sum_tilde_pym2(j)  ;
%         end
%             
%             
            
       
       
        %temp3 is  equivalent to  (y_k -x_j)^2
        temp1 = repmat(datmat(:,m), 1, num_particles)';
        temp2 = repmat(snapshots(:, snaptime_now), 1,N); 
        temp3 = temp1 - temp2;
        %Cap maximum so that we won't divide by zero later.
        temp3 = min(temp3.^2/(2*sigW^2), 745);
        tilde_pym =exp(-temp3);       
        
        %temp3x= temp3.^2/ (2*sigW^2); 
        %centralizer = max(temp3x(:));
        %temp3x = temp3x - centralizer; 
        
        %tilde_pymx = exp(-temp3x); 
        %sum_tilde_pymx = repmat(sum(tilde_pymx,1), num_particles,1);
        %tilde_pymx = tilde_pymx./sum_tilde_pymx;

        
        
        sum_tilde_pym = repmat(sum(tilde_pym,1), num_particles,1);
        tilde_pym = tilde_pym./sum_tilde_pym;
                
        
        deriv_logtemp = permute(repmat(deriv_loglike, 1, 1,num_particles)... 
                        ,[3,1,2]);
        tilde_pytemp = repmat(tilde_pym,1,1,5);
        
        deriv_pym  = deriv_logtemp .* tilde_pytemp; 
   
        deriv_py = mean(deriv_pym,2);
        tilde_py = mean(tilde_pym,2);
        tilde_pys(:,snaptime_now) = tilde_py;
        
        tilde_pytemp2 = repmat(tilde_py,1,1,5);
        deriv_a_frames(:, snaptime_now) = reshape(...
            mean(deriv_py./tilde_pytemp2,1),5,1);
        
        snaptime_now = snaptime_now+1;
      %  toc
    end
    
  end
  
  deriv_a = mean(deriv_a_frames(:,2:end), 2);
  close(h);
%   plots = false;
%   if(plots == true)
%   figure;
%    for(m = 1:N)
%        plot(time, datmat(m,:), 'r'); hold on;
%    end 
%       
%    for(m = 1:N)
%        plot(time, ydat(m,:), 'b'); hold on;
%    end 
%    hold off;      
%       
%   end 


end 

