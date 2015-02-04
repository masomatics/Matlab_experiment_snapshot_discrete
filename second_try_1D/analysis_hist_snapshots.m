%%
%Preamble
D = genpath('../snap1_26_goodfiles');
addpath(D) 
analysis_load_data

%%
%Calling the data, make the histogram structure
[num_particles, num_all_slices] = size(datx)
num_bins = 10; 
field1 = 'weights' ; value1 = zeros(num_bins, num_slices);
field2 = 'values' ; value2 = zeros(num_bins, num_slices);
compress_all_snapshot_hist = struct(field1, value1, field2, value2);

%%
%Convert the data to hisograms
for(k = 1: num_all_slices) 
[num_el, center]= hist(datx(:, k), num_bins);
compress_all_snapshot_hist.weights(:,k) = num_el / sum(num_el);
compress_all_snapshot_hist.values(:,k) = center;
end 
compress_snapshot_hist = ...
    struct(field1, compress_all_snapshot_hist.weights(:,timesample),...
    field2,compress_all_snapshot_hist.values(:,timesample))
%%
%compare the histogram against the real pdf on the state space. 
figure(3);
for(frame = 1: num_slices)
    subplot(2, ceil(num_slices/2), frame); 
    
    %histogram
    plot(compress_snapshot_hist.values(:,frame), compress_snapshot_hist.weights(:,frame)) 
    fig = gca;
    h = findobj(fig,'Type','patch');
    set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
    hold on; 
    scale_h = max(compress_snapshot_hist.weights(:,frame)) ;
    
    %pdf of normal mixture    
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
    scale_pdf = max(pdf)
    plot(xpts, pdf*scale_h/scale_pdf, 'LineWidth', 3)    
    title(['t=', num2str(timesample(frame))])
end 
hold off;%end of histogram plot sequence

%% 
%compare the moments
figure(4);
for(frame = 1 : num_slices)
    m1_data_approx_h(frame) = sum(compress_snapshot_hist.weights(:,frame).*...
        compress_snapshot_hist.values(:,frame).^1);
    m2_data_approx_h(frame) = sum(compress_snapshot_hist.weights(:,frame).*...
        compress_snapshot_hist.values(:,frame).^2);
    m3_data_approx_h(frame) = sum(compress_snapshot_hist.weights(:,frame).*...
        compress_snapshot_hist.values(:,frame).^3);
    m4_data_approx_h(frame) = sum(compress_snapshot_hist.weights(:,frame).*...
        compress_snapshot_hist.values(:,frame).^4);
    m5_data_approx_h(frame) = sum(compress_snapshot_hist.weights(:,frame).*...
        compress_snapshot_hist.values(:,frame).^5);
    m6_data_approx_h(frame) = sum(compress_snapshot_hist.weights(:,frame).*...
        compress_snapshot_hist.values(:,frame).^6);
end 

for(frame = 1 : num_slices)
    subplot(2,3,1)
    plot(timesample, [m1_data_true;m1_data_approx;m1_data_approx_h])
    legend('truth', 'GMM', 'hist')
    title('mean')
%
    subplot(2,3,2)
    plot(timesample, [m2_data_true;m2_data_approx;m2_data_approx_h])
    legend('truth', 'GMM', 'hist')
    title('moment2')
%
    subplot(2,3,3)
    plot(timesample, [m3_data_true;m3_data_approx;m3_data_approx_h]) 
    legend('truth', 'GMM', 'hist')    
    title('moment3')
%
    subplot(2,3,4)
    plot(timesample, [m4_data_true;m4_data_approx;m4_data_approx_h])
    legend('truth', 'GMM', 'hist')
    title('moment4')
%
    subplot(2,3,5)
    plot(timesample, [m5_data_true;m5_data_approx;m5_data_approx_h])  
    legend('truth', 'GMM', 'hist')
    title('moment5')
%
    subplot(2,3,6)
    plot(timesample, [m6_data_true;m6_data_approx;m6_data_approx_h]) 
    legend('truth', 'GMM', 'hist')    
    title('moment6')
    
end 
hold off;
%end of comparison plot sequence