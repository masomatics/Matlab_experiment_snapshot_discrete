% analysis_hist_generation  
%       Creates histogram based compression of the snapshots from the
%       snapshots data.
%
% Input: snapshots[num_particles][num_timesamples]
% Output: structure with attributes weights (histogram weights)
%         and values (histogram supports)             
function compress_snapshot_hist = analysis_hist_generation(snapshots, num_bins)
        %%%Calling the data, make the histogram structure        
        [num_particles, num_slices] = size(snapshots);
        num_bins = 10; 
        field1 = 'weights' ; value1 = zeros(num_bins, num_slices);
        field2 = 'values' ; value2 = zeros(num_bins, num_slices);
        compress_all_snapshot_hist = struct(field1, value1, field2, value2);

        %%%Convert the data to hisograms
        
        for(k = 1: num_slices) 
        [num_el, center]= hist(snapshots(:, k), num_bins);
        compress_all_snapshot_hist.weights(:,k) = num_el / sum(num_el);
        compress_all_snapshot_hist.values(:,k) = center;
        end 
        compress_snapshot_hist = ...
            struct(field1, compress_all_snapshot_hist.weights,...
            field2,compress_all_snapshot_hist.values)

end 