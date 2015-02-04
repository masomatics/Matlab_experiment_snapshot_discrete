function[xdat_test] = debug_GMM_generate(GMM,num_sample)

[numcomp, numpar]= size(GMM);

cumulative = zeros(1,numcomp-1); 
cumulative(1) = GMM(1,3); 
for(m = 2 : numcomp-1)
    cumulative(m) = cumulative(m-1) + GMM(m-1,3); 
end 



unifseed = rand(1, num_sample);
randseed = randn(1, num_sample);
xdat_test = NaN(1,num_sample);
    for(k = 1 : num_sample)

        component  = 1; 
        for(m = 1 :numcomp-1)
            if(unifseed(k) > cumulative(m))
            component = component  + 1;
            end
        end 
        xdat_test(k) = randseed(k) *  sqrt(GMM(component,2)) +  GMM(component,1);


    end 

end 
