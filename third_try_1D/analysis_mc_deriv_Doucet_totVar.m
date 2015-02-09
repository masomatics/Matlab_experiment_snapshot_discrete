%analysis_mc_deriv_Doucet_totalV
%
%This computes the derivative of the objective function  (Total Variation
%distance)
% |p( dot, theta)- q( dot )|
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
%   rnsource_sys    :std normal randoms from which the path of x will be generated.
%   compress_snap_vals: support of the approx. empirical distribution  
%   compress_snap_wgts: weight on the approx. empirical distribution
%   snaptime    :times at which the snapshots were taken. The length is
%               num_frames.
%   N           :number of Xpaths to be sampled
%
% Outputs: 
%
%   derivative  :Derivative of the energy wrt parameters  T by num_parameters. 
%   energy      :E_q[log p]
%   p_ymk_all   :the approx of p in 2.(a).(iv) of pseudocode. num_particles
%                by num_frames.
%   


function [derivative, energy, p_ymk_all] = analysis_mc_deriv_Doucet_totVar(initmean, ...
theta, sigV, sigW, T, rnsource_sys, compress_snap_vals, compress_snap_wgts, snaptime, N)

%%%preset variables%%%
[num_particles, num_slices] = size(compress_snap_vals); 
 num_parameters = length(theta);   
include = [];
%%%Target variables%%%
    derivative= zeros(num_parameters,1); %Gradient of the parameter
    tilde_p_ymk = zeros(num_particles);
    p_ymkj = zeros(N, num_particles);
    dp_jr = zeros(N, num_particles);
    dEP_kr = zeros(N, num_particles); 
    
%%%Temporary variables%%%
    datmat = NaN(1,N);
    ftheta = NaN(5,N);
    ycopy = zeros(num_particles, N); %used in comparing all y to x.
    xcopy = zeros(num_particles, N); %used in comparing all y to x.
    scale_p_ymkj = zeros(1,N);
    slace_p_ymkj_mat = zeros(num_particles, N); 
    derivative_alltime = zeros(num_parameters, num_slices); %derivative at each time. 
%%%Initialization%%%
    v = sigV*rnsource_sys;  %system noise
    initx = sqrt(5)* v(1,:)  + initmean; %common initial distrbn
    datmat = initx;  %initial value of the system 
    deriv_loglike = zeros(num_parameters,N);
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
            p_ymkj= exp(-min((ycopy - xcopy).^2/sigW^2,300)); %(5.45)
            scale_p_ymkj = 1./sum(p_ymkj,1);
            scale_p_ymkj_mat = repmat(scale_p_ymkj, [num_particles,1]); 
            p_ymkj = p_ymkj .* scale_p_ymkj_mat;
            
            tilde_p_ymk = mean(p_ymkj,2);
            dEP_kr = 1/N * p_ymkj * deriv_loglike';
            
            derivative_alltime(:,snaptime_now) = ...
                sign(tilde_p_ymk -compress_snap_wgts(:,snaptime_now))' * ...
                (dEP_kr);            
            
            %DEBUG : Make sure that the following two plots are the SAME
            %when using the same noise and same theta. 
            %plot(compress_snap_vals(:,snaptime_now),compress_snap_wgts(:,snaptime_now)*N)
            %hold on;
            %[hoge, piyo] =hist(datmat);
            %plot(piyo, hoge)
            %hold off;
            p_ymk_all(:,snaptime_now) = tilde_p_ymk;
            %DEBUG:Compute energy ; Debug
            %energy(snaptime_now) = sum(log(p_ymk/sum(p_ymk)).* ...
            %    compress_snap_wgts(:,snaptime_now));
            energy(snaptime_now) = sum((tilde_p_ymk -compress_snap_wgts(:,snaptime_now)).^2);
                   if(max(p_ymkj(:))  > min(p_ymkj(:)) )
                        %In the case of using empirical distribution,
                        %the true distribution can be soo far from the
                        %simulated distribution that the energy can be
                        %wrongfully too high (all uniform) This part is
                        %very hand-wavy.... I am ridding of such frames
                        %from the consideration.                         
                        include = [include, snaptime_now];
                   end          
            snaptime_now = snaptime_now +1;
        end 



    end 
    derivative = sum(derivative_alltime(:,(include-1)),2);
end


function datmat_new = internal_update(datmat, theta, noise)
        datmat_new = theta(1) + theta(2).*datmat + ...
            theta(3).*datmat./ (1 + (datmat.^2)) + ...
            theta(4).*cos(theta(5).*datmat) + noise ;
end 
