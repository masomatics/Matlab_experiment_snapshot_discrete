% Generate a size N batch of simulation paths 
% This for the discrete process used in the paper of Doucet. 
%
%
% x(m) = theta(1) + theta(2)*x(m-1) + ...
%         theta(3)*x(m-1)/ (1 + x(m-1)^2) + ...
%         theta(4)*cos(theta(5)*x(m-1)) + N(0, sigV^2)
% y(m) = f(x(m) , b)  + N(0, sigW^2)
%
% Input 
%   double[dimen] initmean
%   int  T >0
%   int N  batchsize
%   Covariance Matrix sigV
%   Covariance Matrix sigW
%   int initmean  initial point mean
%   double[N][num_timepts] rnsource_sys  
%   double[N][num_timepts] rnsource_obs  
%   double[num_components, 3, proportion] compress_snapshots 
%
% Output 
%   double[num_slices] sensitivity
%



function [val_sensitivity, energies] = analysis_sensitivity(compress_data, snaptime,...
    initmean,theta, T, sigV,  rnsource_sys, N)

  %default = true;
  default = false;
  if(default == true)
    sigV = sqrt(10);
  end
  
  num_parameters = length(theta);
  
  [num_components, num_GMMparameters, num_frames] = size(compress_data);
  tnow = 1;
  snapframe_now = 2;
  temp_xdat = NaN(N,T);
  xdat = NaN(N,1);
  deriv_log_p = zeros(N,5);
  delta_deriv_log_p = zeros(N,5);
  val_sensitivity = zeros(num_frames,num_parameters);  
  objectiv = zeros(N,1);
  v = sigV*rnsource_sys;
  initx = sqrt(5)* randn(N,1)  + initmean;
  energy = 0;
  energies = zeros(1, num_frames);
  %initx = initmean;    
  xdat = initx;
  
  
  for(k = 2:T)
    %xdat - xdat_hat = v(:,k)
    %-(xdat - xdat_hat)* deriv(xdat_hat) is the answer
%  
%       deriv_log_p2(:,1) = deriv_log_p(:,1) - v(:,k).* 1./sigV;
%       deriv_log_p2(:,2) = deriv_log_p(:,2) - v(:,k).* xdat ./sigV;
%       deriv_log_p2(:,3) = deriv_log_p(:,3) - v(:,k).* xdat ./ (1 + (xdat).^2) ./sigV;
%       deriv_log_p2(:,4) = deriv_log_p(:,4) - v(:,k).* cos(theta(5).*xdat) ./sigV;
%       deriv_log_p2(:,5) = deriv_log_p(:,5) - v(:,k).* (-theta(4).*xdat ...
%           .*sin(theta(5).*xdat))./sigV;    

    delta_deriv_log_p = -repmat(v(:,k), 1,5) .*... 
    [ ones(N,1), xdat, xdat ./ (1 + (xdat).^2) ,...
     cos(theta(5).*xdat), (-theta(4).*xdat.*sin(theta(5).*xdat))]/sigV;

    deriv_log_p = deriv_log_p + delta_deriv_log_p;

    xdat = theta(1) + theta(2)*xdat+ ...
         theta(3)*xdat./ (1 + xdat.^2) + ...
         theta(4)*cos(theta(5)*xdat) + v(:,k) ;   
    
    if(k == snaptime(snapframe_now))
        % 
        for(m = 1 : num_components)
            objectiv = objectiv + compress_data(m, 3, snapframe_now) * ...
                exp(-(xdat - compress_data(m, 1, snapframe_now) ).^2/ ...
                (2*compress_data(m, 2, snapframe_now))...
                )./sqrt(compress_data(m, 2, snapframe_now)); 
        end
       val_sensitivity(snapframe_now, :) = mean((repmat(log(objectiv), 1,num_parameters) .* ...
            deriv_log_p),1);
        energies(snapframe_now) = mean(log(objectiv));       
        objectiv = zeros(N,1);

                            num_slices = length(snaptime);
                            frame = snapframe_now;
                            subplot(2, ceil(num_slices/2), frame)
                            num_breaks = 50;
                            hist(xdat,num_breaks);
                            fig = gca;
                            h = findobj(fig,'Type','patch');
                            set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
                            hold on; 
                            limit1 = xlim(fig);

                            num_xbreaks = 1000;
                            xpts = limit1(1) + (limit1(2) - limit1(1)) * (1:num_xbreaks)/ num_xbreaks;

                            pdf = zeros(1, length(xpts));
                            for(compo = 1: num_components)
                                if(compress_data(compo, 3, frame) > 0)
                                    pdf = pdf + compress_data(compo, 3, frame) * ...
                                        1/sqrt(2*pi*(compress_data(compo, 2, frame)))*...
                                        exp(-(xpts - compress_data(compo, 1, frame)).^2./ ...
                                        (2*compress_data(compo, 2, frame))) ;
                                end
                            end 
                            plot(xpts, pdf*N, 'LineWidth', 3)

                            title(['t=', num2str(snaptime(snapframe_now))])        
        
        snapframe_now = min(snapframe_now  + 1, length(snaptime));

    end

        
 
        
    
  end
end 

