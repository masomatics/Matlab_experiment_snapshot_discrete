%% Reload data and fit Gaussian 
close all;

theta_a = [0, 1/2, 25, 6, 0.2];
theta_b = 20;

sigV = sqrt(10);
sigW = sqrt(1); 


init = -1.5


datx= load('/Users/Markov/Desktop/SnapShot_Diffusion2/October_Offensive/SnapShot_Matlabcode_beta4/xdata6000particles.mat');
daty= load('/Users/Markov/Desktop/SnapShot_Diffusion2/October_Offensive/SnapShot_Matlabcode_beta4/ydata6000particles.mat');


datx= datx.xdatm;
daty= daty.ydatm;

timesample = [1, 5:5:40];

snapshots = daty(:,timesample);
dimension = size(snapshots) ;
num_pts = dimension(1);
num_slices = dimension(2);
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
        AIC= GMModel.AIC;
        new_num_components = max_num_components; 
        while(gosign == false)
           num_try = num_try + 1; 
           display([num2str(num_try), 'st try, ', num2str(new_num_components), ...
               'components' ]);
           newGMModel = gmdistribution.fit(snapshots(:,k), new_num_components);
           gosign1= newGMModel.Converged;
           
           gosign = gosign1;
           display([AIC, newAIC]);

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
for(frame = 1: num_slices)
    subplot(2, ceil(num_slices/2), frame)
    num_breaks = 50;
    hist(snapshots(:,frame),num_breaks);
    fig = gca;
    h = findobj(fig,'Type','patch');
    set(h,'FaceColor','r','EdgeColor','w','facealpha',0.5)
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
    plot(xpts, pdf*num_pts, 'LineWidth', 3)

end
filename = ['GMMfit_graph_' , num2str(max_num_components), '_components.jpg' ]
filecommit =['GMMfit_graph_' , num2str(max_num_components), '_components_commit.jpg' ]

saveas(hoge,filename);
%saveas(hoge,filecommit);

% 
% close all;
% numsample = 6000
% varian = 2;
% xx = randn(numsample,1)* sqrt(varian);
% hist(xx, num_breaks)
% hold on;
% figfig = gca; 
% limitx = xlim(figfig);
% xxpts = limitx(1) + (limitx(2) - limitx(1)) * (1:num_breaks)/ num_breaks;
% pdfpdf = 1/sqrt(2*pi*varian)*exp(-(xxpts - 0).^2./(2*varian)) ;
% plot(xxpts, pdfpdf*numsample)




