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


function [derivative, derivative_alltime, statmean] = debug_analysis_mc_deriv_simplestat(initmean, ...
theta, sigV, T, rnsource_sys, snaptime, N)

%%%Preset Variables %%%
    num_slices = length(snaptime);
    num_parameter = length(theta);
%%%Target variables%%%
    derivative= zeros(num_parameter,1); %Gradient of the parameter 
    statmean = zeros(1, num_slices);
%%%Temporary variables%%%
    datmat = NaN(1,N);
    ftheta = NaN(1,N,num_parameter);
    xcopy = zeros(N,num_parameter); %used in comparing all y to x.
    derivative_alltime = zeros(5, num_slices); %derivative at each time. 
%%%Initialization%%%
    v = sigV*rnsource_sys;  %system noise
    initx = sqrt(5)* v(1,:) + initmean; %common initial distrbn
    datmat = initx;  %initial value of the system 
    deriv_loglike = zeros(1,N,num_parameter);
    snaptime_now = 2;     
%%%Monte Carlo Simulation%%% 
    for(time = 2:T)
        %compute deriv f% 
        ftheta(1,:,1) = ones(1,N);
        ftheta(1,:,2) = datmat;
        ftheta(1,:,3) = datmat./ (1 + (datmat).^2);
        ftheta(1,:,4) = cos(theta(5).*datmat);
        ftheta(1,:,5) = theta(4).*datmat.*(-sin(theta(5).*datmat));
        %Update log likelihood : kth component is partial k, (5.44)%
        deriv_loglike = deriv_loglike + ...
            repmat(v(time,:),1,1,num_parameter).*ftheta / sigV^2 ;         
        %Update X_t%
        datmat = internal_update(datmat, theta, v(time,:));         
        %%DEBUG
        %daty = datmat + w(time,:);        
        if(time == snaptime(snaptime_now))
            %
            xcopy = repmat(datmat,num_parameter,1);
            %
            dp_xjmr= xcopy .* squeeze(deriv_loglike)'; %(5.46)
            dp_xmr = squeeze(mean(dp_xjmr,2));  %(5.47) size (num_particles * 5) 
            % (5.41)      
            derivative_alltime(:,snaptime_now) = dp_xmr;
            statmean(snaptime_now) = mean(datmat);
            snaptime_now = snaptime_now + 1;
        end 
    end 
    derivative = sum(derivative_alltime,2);
end


function datmat_new = internal_update(datmat, theta, noise)
        datmat_new = theta(1) + theta(2).*datmat + ...
            theta(3).*datmat./ (1 + (datmat.^2)) + ...
            theta(4).*cos(theta(5).*datmat) + noise ;
end 
