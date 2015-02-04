num_sample  = 10000

GMMtruth = compress_snapshots(:,:,3);
GMMfake1 = compress_snapshots(:,:,1);
% [numcomp, numpar]= size(compress_snapshots(:,:,2));
% cumulative = zeros(1,numcomp-1); 
% cumulative(1) = GMM(1,3); 
% for(m = 2 : numcomp-1)
%     cumulative(m) = cumulative(m-1) + GMM(m-1,3); 
% end 
% 
% 
% 
% unifseed = rand(1, num_sample);
% randseed = randn(1, num_sample);
% xdat_test = NaN(1,num_sample);
% for(k = 1 : num_sample)
%     
%     component  = 1; 
%     for(m = 1 :numcomp-1)
%         if(unifseed(k) > cumulative(m))
%         component = component  + 1;
%         end
%     end 
%     xdat_test(k) = randseed(k) *  sqrt(GMM(component,2)) +  GMM(component,1);
%     
%     
% end 
xdat_truth = debug_GMM_generate(GMMtruth ,num_sample);
xdat_fake0 = randn(1,num_sample)*1.5;
xdat_fake1 = debug_GMM_generate(GMMfake1 ,num_sample);




objectiv_truth = 0; 
 for(m = 1 : num_components)
            objectiv_truth = objectiv_truth + GMMtruth(m, 3) * ...
                exp(-(xdat_truth - GMMtruth(m, 1) ).^2/ ...
                (2*GMMtruth(m, 2))...
                )./sqrt(GMMtruth(m, 2)); 
 end
 mean(log(objectiv_truth))
 
 objectiv_fake0= 0;
for(m = 1 : num_components)
            objectiv_fake0 = objectiv_fake0 + GMMtruth(m, 3) * ...
                exp(-(xdat_fake0 - GMMtruth(m, 1) ).^2/ ...
                (2*GMMtruth(m, 2))...
                )./sqrt(GMMtruth(m, 2)); 
 end
 mean(log(objectiv_fake0))
 
 
 objectiv_fake1= 0;
for(m = 1 : num_components)
            objectiv_fake1 = objectiv_fake1 + GMMtruth(m, 3) * ...
                exp(-(xdat_fake1 - GMMtruth(m, 1) ).^2/ ...
                (2*GMMtruth(m, 2))...
                )./sqrt(GMMtruth(m, 2)); 
 end
 mean(log(objectiv_fake1))
 