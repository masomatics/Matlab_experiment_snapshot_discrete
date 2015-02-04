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
%   Covariance Matrix sigV
%   Covariance Matrix sigW
%   double[N][num_timepts] rnsource_sys  
%   double[N][num_timepts] rnsource_obs  
%   int N number of experiments
%
% Output 
%   xdat
%   ydat 
%   tdat



function [xdat, ydat] = analysis_data_generation_N(initmean,theta, T, sigV, sigW,  rnsource_sys, rnsource_obs, N)

  %default = true;
      b = 20;

  default = false;
  if(default == true)
    %N = 100;
%    a = [1/2, 25, 8, 1.2];
    theta = [0, 1/2, 25, 6, 0.2];
    b = 20;
%    b =10;   
    sigV = sqrt(10);
    sigW = sqrt(2);
%   sigW = sqrt(5);
  end
  
  parameters = {theta, b, sigV, sigW};
  display(parameters);
  
  
  time =  1:T  ;
  xdat = NaN(N,T);
  ydat = NaN(N,T);
  datmat= NaN(2,N,T);
  
  
  v = sigV*rnsource_sys;
  w = sigW*rnsource_obs;
 initx = sqrt(5)* v(:,1)  + initmean;
 %initx = initmean;
  
  %initx = sqrt(120);
 % inity = 1/b*initx.^(2) + w(:,1); 
  inity = initx + w(:,1);
  
  xdat(:,1) = initx;
  ydat(:,1) = inity;
  
  
  for(k = 2:T)
     xdat(:,k) = theta(1)+ theta(2).*xdat(:,k-1) + theta(3).*xdat(:,k-1)./ (1 + (xdat(:,k-1)).^2) + ...
         theta(4).*cos(theta(5).*xdat(:,k-1)) + v(:,k);
     %ydat(:,k) = 1/b .* xdat(:,k).^2 + w(:,k);  
     ydat(:,k) = xdat(:,k) + w(:,k);  

  end
  
  plots = false;
  if(plots == true)
  figure;
   for(m = 1:N)
       plot(time, xdat(m,:), 'r'); hold on;
   end 
      
   for(m = 1:N)
       plot(time, ydat(m,:), 'b'); hold on;
   end 
   hold off;      
      
  end 


end 

