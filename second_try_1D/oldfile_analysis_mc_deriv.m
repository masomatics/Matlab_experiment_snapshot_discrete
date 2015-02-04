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


function [derivative] = analysis_mc_deriv(initmean, ...
theta, sigV, sigW, T, rnsource, compress_snap_vals, compress_snap_wgts, snaptime, N)

%%%preset variables%%%
[num_particles, num_slices] = size(compress_snap_vals); 
 num_parameters = 5;   

%%%Target variables%%%
    derivative= zeros(5,1); %Gradient of the parameter
    p_ymk = zeros(num_particles,1);   %2.(c).ii, fixed k.
    dp_ymkr = zeros(num_particles,5); % 2.(c).ii fixed k.
    
%%%Temporary variables%%%
    datmat = NaN(1,N);
    ftheta = NaN(1,N,5);
    p_ymkj = zeros(num_particles, N);    %2.(c).i fixed k.
    deriv_p_ymkj = zeros(num_particles, N,5); %2.(c).i  fixed k.    
    ycopy = zeros(num_particles, N); %used in comparing all y to x.
    xcopy = zeros(num_particles, N); %used in comparing all y to x.
    xypair= zeros(num_particles, N); %N(y|x,sigV)
    deriv_loglike_copy= zeros(num_particles,N,5);
    p_ymkj_copy = zeros(num_particles, N,5);
    p_ymk_copy = zeros(num_particles,5);
    dp_ymkjr = zeros(num_particles, N, 5);
    derivative_alltime = zeros(5, num_slices); %derivative at each time. 
%%%Initialization%%%
    v = sigV*rnsource;  %system noise
    initx = sqrt(5)*v(1,:)  + initmean; %common initial distrbn
    datmat = initx;  %initial value of the system 
    deriv_loglike = zeros(1,N,5);
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
            repmat(v(time,:),1,1,5).*ftheta / sigV^2 ; 

        %Update X_t%
        datmat = theta(1) + theta(2).*datmat + ...
            theta(3).*datmat./ (1 + (datmat.^2)) + ...
            theta(4).*cos(theta(5).*datmat) + v(time,:) ;
        
        if(time == snaptime(snaptime_now))
            ycopy = repmat(compress_snap_vals(:,snaptime_now), 1,N);
            xcopy = repmat(datmat,num_particles,1);
            p_ymkj= exp(-(ycopy - xcopy).^2/sigW^2); %(5.45)
            p_ymkj_copy = repmat(p_ymkj, [1,1,5]);
            deriv_loglike_copy = repmat(deriv_loglike,num_particles,1,1);
            dp_ymkjr= p_ymkj_copy.* deriv_loglike_copy; %(5.46)
            p_ymk = mean(p_ymkj,2);  %(5.47) size (num_particles,1)
            dp_ymkr = squeeze(mean(dp_ymkjr,2));  %(5.47) size (num_particles * 5) 
            p_ymk_copy = repmat(...
            compress_snap_wgts(:,snaptime_now).*...
               (1./p_ymk - 1/sum(p_ymk.*compress_snap_wgts(:,snaptime_now))),1,5);
             % (5.41)      
            derivative_alltime(:,snaptime_now) = sum(dp_ymkr.*p_ymk_copy...
                , 1);             
            snaptime_now = snaptime_now +1;
        end 
    end 
    derivative = sum(derivative_alltime,2);
end


