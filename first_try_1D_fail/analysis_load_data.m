close all;
theta = [0, 1/2, 25, 6, 0.2];
%theta = [-1.5,0.1,0,0,0]; 

num_timepts = 40;
sys_noise = randn(1, num_timepts);
obs_noise = randn(1, num_timepts);
sigV = sqrt(10);
sigW = sqrt(2); 
init = -1.5

N = 6000;

sys_noise = randn(N,num_timepts);
obs_noise = randn(N,num_timepts);
%% 

%[datx1, daty1]= oldfile_data_generation_N2(init, theta, num_timepts, sigV, sigW,  sys_noise, obs_noise , N )

[datx, daty]= analysis_data_generation_N(init, theta, num_timepts, sigV, sigW,  sys_noise, obs_noise , N );

%%

timesample = [1, 5:5:40];
%timesample = [1,10];

snapshots = daty(:,timesample);



%%

max_num_components = 5;
AIC_epsilon = 0.0001
AIC_option = false
%(1) Mean 
%(2) Sigma
%(3) Weight
compress_snapshots = zeros(max_num_components, 3, length(timesample));


%Fit Mixture of Gaussians AIC choice mode
if(AIC_option == true)
    for(k = 1: length(timesample))
        display(['fitting the snap number ' , num2str(k)]) 
        gosign = false;
        num_try = 0;
        num_components = 1;
        GMModel = gmdistribution.fit(snapshots(:,k), 1);
        AIC= GMModel.AIC;
        new_num_components = num_components +1; 
        while(gosign == false)
           num_try = num_try + 1; 
           display([num2str(num_try), 'st try, ', num2str(new_num_components), ...
               'components' ]);
           newGMModel = gmdistribution.fit(snapshots(:,k), new_num_components);
           newAIC = newGMModel.AIC;
           gosign1= newGMModel.Converged;
           gosign2= (newAIC > AIC);
           if(gosign1*(1- gosign2))
               num_components = new_num_components;
               oldAIC = AIC;
               AIC = newAIC;
               GMModel = newGMModel;
               new_num_components = new_num_components+1;
           end                

           gosign = gosign1*gosign2;
           display([AIC, newAIC]);

        end 
        display(['frame', num2str(k), ' success with ', ...
            num2str(num_components), ' components']); 
        compress_snapshots(1:num_components,1,k) =  GMModel.mu ;
        compress_snapshots(1:num_components,2,k) =  GMModel.Sigma(:) ;
        compress_snapshots(1:num_components,3,k) =  GMModel.PComponents;
    end
end

if(AIC_option == false)
%%Fit Mixture of Gaussians: Max mixture mode:
    for(k = 1: length(timesample))
        display(['fitting the snap number ' , num2str(k)]) 
        gosign = false;
        num_try = 0;
        num_components = max_num_components;
        GMModel = gmdistribution.fit(snapshots(:,k), 1);
        new_num_components = max_num_components; 
        while(gosign == false)
           num_try = num_try + 1; 
           display([num2str(num_try), 'st try, ', num2str(new_num_components), ...
               'components' ]);
           newGMModel = gmdistribution.fit(snapshots(:,k), new_num_components);
           gosign1= newGMModel.Converged;
           
           gosign = gosign1;
        end 
        GMModel = newGMModel
        display(['frame', num2str(k), ' success with ', ...
            num2str(num_components), ' components']); 
        compress_snapshots(1:num_components,1,k) =  GMModel.mu ;
        compress_snapshots(1:num_components,2,k) =  GMModel.Sigma(:) ;
        compress_snapshots(1:num_components,3,k) =  GMModel.PComponents;
    end 
end


%%
hoge = figure(1);
%Number of breaks will change the height dependency. Graphs might differ in
%heights. 
num_slices = length(timesample)
for(frame = 1: num_slices)
    subplot(2, ceil(num_slices/2), frame)
    num_breaks = 50;
    hist(snapshots(:,frame),num_breaks);
    fig = gca;
    h = findobj(fig,'Type','patch');
    set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
    hold on; 
    limit1 = xlim(fig);

    num_xbreaks = 1000;
    xpts = limit1(1) + (limit1(2) - limit1(1)) * (1:num_xbreaks)/ num_xbreaks;

    pdf = zeros(1, length(xpts));
    for(compo = 1: max_num_components)
        if(compress_snapshots(compo, 3, frame) > 0)
            pdf = pdf + compress_snapshots(compo, 3, frame) * ...
                1/sqrt(2*pi*(compress_snapshots(compo, 2, frame)))*...
                exp(-(xpts - compress_snapshots(compo, 1, frame)).^2./ ...
                (2*compress_snapshots(compo, 2, frame))) ;
        end
    end 
    plot(xpts, pdf*N, 'LineWidth', 3)
    
    title(['t=', num2str(timesample(frame))])
end




%% Compare the mean, variance and higher moments 

%Mean, variance

m1_data_true =  mean(snapshots, 1);
var_data_true= var(snapshots, 1);
m3_data_true = mean(snapshots.^3, 1);
m4_data_true = mean(snapshots.^4, 1);
m5_data_true = mean(snapshots.^5, 1);
m6_data_true = mean(snapshots.^6, 1);

summary_truth = [m1_data_true; var_data_true; m3_data_true ; ...
    m4_data_true; m5_data_true; m6_data_true] ;



m1_data_approx =  zeros(1,9);
m2_data_approx = zeros(1,9);
var_data_approx = zeros(1,9);
m3_data_approx =zeros(1,9);
m4_data_approx =zeros(1,9);
m5_data_approx =zeros(1,9);
m6_data_approx =zeros(1,9);


for(frame = 1 : num_slices)
    for(compo = 1 : num_components)
        m1_data_approx(frame) = m1_data_approx(frame)+ ...
            compress_snapshots(compo,1,frame) * ...
            compress_snapshots(compo,3,frame);
        m2_data_approx(frame) = m2_data_approx(frame)+ ...
            (compress_snapshots(compo,2,frame) + compress_snapshots(compo,1,frame)^2 )  * ...
            compress_snapshots(compo,3,frame);
        m3_data_approx(frame) = m3_data_approx(frame)+ ...
            (compress_snapshots(compo,1,frame)^3 + ....
            3*compress_snapshots(compo,1,frame)*compress_snapshots(compo,2,frame)) *...
            compress_snapshots(compo,3,frame);
        m4_data_approx(frame) = m4_data_approx(frame)+ ...
            (compress_snapshots(compo,1,frame)^4 + ...
            6*compress_snapshots(compo,1,frame)^2*compress_snapshots(compo,2,frame)+ ...
            3*compress_snapshots(compo,2,frame)^2) * ...
            compress_snapshots(compo,3,frame);
        m5_data_approx(frame) = m5_data_approx(frame)+ ...
            (compress_snapshots(compo,1,frame)^5 + ...
            10*compress_snapshots(compo,1,frame)^3*compress_snapshots(compo,2,frame)+...
            15*compress_snapshots(compo,1,frame)*compress_snapshots(compo,2,frame)^2) * ...
            compress_snapshots(compo,3,frame);
        m6_data_approx(frame) = m6_data_approx(frame)+ ...
            (compress_snapshots(compo,1,frame)^6 + ...
            10*compress_snapshots(compo,1,frame)^4*compress_snapshots(compo,2,frame)+...
            45*compress_snapshots(compo,1,frame)^2*compress_snapshots(compo,2,frame)^2 + ...
            15*compress_snapshots(compo,2,frame)^3) * ...
            compress_snapshots(compo,3,frame);            
        
    end    
    var_data_approx(frame) = m2_data_approx(frame) - m1_data_approx(frame)^2;
end

figure(2)

subplot(2,3,1)
plot(timesample, m1_data_true)
hold on;
plot(timesample, m1_data_approx)
title('mean')


subplot(2,3,2)

plot(timesample, var_data_true)
hold on;
plot(timesample, var_data_approx)
title('var')


subplot(2,3,3)

plot(timesample, m3_data_true)
hold on;
plot(timesample, m3_data_approx)
title('moment3')



subplot(2,3,4)

plot(timesample, m4_data_true)
hold on;
plot(timesample, m4_data_approx)
title('moment4')


subplot(2,3,5)

plot(timesample, m5_data_true)
hold on;
plot(timesample, m5_data_approx)
title('moment5')


subplot(2,3,6)

plot(timesample, m6_data_true)
hold on;
plot(timesample, m6_data_approx)
title('moment6')

%% Compare even the higher moments







