function half_peak_distance = Find_Half_Peak(dist, amp)

%%find the half peak points:
%% to smooth first
filter_size = 10;
for i = 1: length(dist)
    amp_filtered(i) = median(amp(max(1, (i-ceil(filter_size/2))):min((i+ceil(filter_size/2)),length(amp))));
end
% plot(dist, amp, 'k-'); hold on
% plot(dist, amp_filtered, 'r*-');

left_index = find(dist <0); %%left to the peak
point_index = find(amp_filtered (left_index) > 0.5); 
left_distance = abs(dist(left_index(point_index(1)))); %% the first point above 50%

right_index = find(dist >0); %%left to the peak
point_index = find(amp(right_index)> 0.5);
right_distance = dist(right_index(point_index(end))); %% the last point above 50%

half_peak_distance = left_distance + right_distance;