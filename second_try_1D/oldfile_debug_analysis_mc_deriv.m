%analysis_mc_deriv
%
%This computes the derivative of the objective function 
%E_q [log_\theta p(y, theta)]
%when q is an approximation of the given snapshot. 
%The algorithm follows 5.6 of my snapnote. 
%
%
%
% The model of interest: 
%
% x(m) = theta(1) + theta(2)*x(m-1) + ...
%         theta(3)*x(m-1)/ (1 + x(m-1)^2) + ...
%         theta(4)*cos(theta(5)*x(m-1)) + N(0, sigV^2)
% y(m) = N(x(m), sigW^2);
% x(1) = N(initmean, sigW^2)
%
% Inputs: 
%   initmean    :As described above
%   theta       :As described above
%   sigV        :System noise.
%   sigW        :Observation noise. 
%   T           :Terminal time. Discrete.
%   rnsource    :std normal randoms from which the path of x will be generated.
%   compress_snap_vals: support of the approx. empirical distribution  
%   compress_snap_wgts: weight on the approx. empirical distribution
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


function [derivative, energy, p_ymk_all] = debug_analysis_mc_deriv(initmean, ...
theta, sigV, sigW, T, rnsource_sys, rnsource_obs, compress_snap_vals, compress_snap_wgts, snaptime, N)

%%%preset variables%%%
[num_particles, num_slices] = size(compress_snap_vals); 
 num_parameters = 5;   

%%%Target variables%%%
    derivative= zeros(num_parameters,1); %Gradient of the parameter
    p_ymk = zeros(num_particles,1);   %2.(c).ii, fixed k.
    dp_ymkr = zeros(num_particles,num_parameters); % 2.(c).ii fixed k.
    
%%%Temporary variables%%%
    datmat = NaN(1,N);
    ftheta = NaN(5,N);
    p_ymkj = zeros(num_particles, N);    %2.(c).i fixed k.
    ycopy = zeros(num_particles, N); %used in comparing all y to x.
    xcopy = zeros(num_particles, N); %used in comparing all y to x.
    derivative_alltime = zeros(num_parameters, num_slices); %derivative at each time. 
%%%Initialization%%%
    v = sigV*rnsource_sys;  %system noise
    w = sigW*rnsource_obs;
    initx = sqrt(5)* v(1,:)  + initmean; %common initial distrbn
    datmat = initx;  %initial value of the system 
    deriv_loglike = zeros(5,N);
    snaptime_now = 2; 
    
%%%DEBUG: Assessment variable%%%
    energy = zeros(1, num_slices);  
    p_ymk_all = zeros(num_particles, num_slices); 
    
%%%Monte Carlo Simulation%%% 
    for(time = 2:T)
    
        %compute deriv f% 
        ftheta(1,:) = ones(1,N);
        ftheta(2,:) = datmat;
        ftheta(3,:) = datmat./ (1 + (datmat).^2);
        ftheta(4,:) = cos(theta(5).*datmat);
        ftheta(5,:) = theta(4).*datmat.*(-sin(theta(5).*datmat));
        %Update log likelihood : kth component is partial k, (5.44)%
        deriv_loglike = deriv_loglike + ...
            repmat(v(time,:),num_parameters,1).*ftheta / sigV^2 ;         
        %Update X_t%
        datmat = internal_update(datmat, theta, v(time,:));         
        %%DEBUG
        %daty = datmat + w(time,:);
        
        if(time == snaptime(snaptime_now))
            ycopy = repmat(compress_snap_vals(:,snaptime_now), 1,N);
            xcopy = repmat(datmat,num_particles,1);
            p_ymkj= exp(-(ycopy - xcopy).^2/sigW^2); %(5.45)
            p_ymk = mean(p_ymkj,2);  %(5.47) size (num_particles,1)
            
            dp_ymkr = 1/N*p_ymkj * deriv_loglike';

            derivative_alltime(:,snaptime_now)=(compress_snap_wgts(:,snaptime_now).*...
               (1./p_ymk - 1/sum(p_ymk.*compress_snap_wgts(:,snaptime_now))))'...
               *dp_ymkr;            
            
            %DEBUG : Make sure that the following two plots are the SAME
            %when using the same noise and same theta. 
            %plot(compress_snap_vals(:,snaptime_now),compress_snap_wgts(:,snaptime_now)*N)
            %hold on;
            %[hoge, piyo] =hist(datmat);
            %plot(piyo, hoge)
            %hold off;
            
            %DEBUG:Compute energy ; Debug
            energy(snaptime_now) = sum(log(p_ymk).* ...
                compress_snap_wgts(:,snaptime_now));
            p_ymk_all(:,snaptime_now) = p_ymk/ sum(p_ymk); 
            snaptime_now = snaptime_now +1;
        end 



    end 
    derivative = sum(derivative_alltime,2);
end


function datmat_new = internal_update(datmat, theta, noise)
        datmat_new = theta(1) + theta(2).*datmat + ...
            theta(3).*datmat./ (1 + (datmat.^2)) + ...
            theta(4).*cos(theta(5).*datmat) + noise ;
end 
