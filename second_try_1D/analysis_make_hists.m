function compress_snapshot_hist = analysis_make_hists...
    (datx, num_bins, num_slices,num_all_slices, timesample)

field1 = 'weights' ; value1 = zeros(num_bins, num_slices);
field2 = 'values' ; value2 = zeros(num_bins, num_slices);
compress_all_snapshot_hist = struct(field1, value1, field2, value2);
for(k = 1: num_all_slices) 
[num_el, center]= hist(datx(:, k), num_bins);
compress_all_snapshot_hist.weights(:,k) = num_el / sum(num_el);
compress_all_snapshot_hist.values(:,k) = center;
end 
compress_snapshot_hist = ...
    struct(field1, compress_all_snapshot_hist.weights(:,timesample),...
    field2,compress_all_snapshot_hist.values(:,timesample))
end