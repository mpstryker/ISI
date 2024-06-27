
function IMAPhase = IMA_draw(line_select)
% Matlab code for map computation
% For ploting phase along aline or contour
% Jianhua Cang, 12-21-04

% important note: x is colomn, and y is row
% Matlab use M(ro, co), i.e., M(y,x)

global mm_per_pixel deg_per_stimcycle filter_size;
global TemplateMethod TemplateValue;
global MapAnsColorMap;
global IMAamp;

switch line_select

case 'contour' % call makecontour.m
    contour_value = str2num(get(findobj('tag','contourvalueeditnew'), 'string'));
    contour_structure = makecontour(contour_value);
    plot_y = contour_structure(1).isoY; 
    plot_x = contour_structure(1).isoX;
    [plot_x, index]= sort(plot_x);
    plot_y = plot_y(index);
    %read in the target map    
    [map_strength, map_phase, file_name] = read_map('Please select the target map', filter_size);
    map_phase = map_phase/360 *deg_per_stimcycle; % deg in real distance
    hndl = figure;
    subplot(2, 1, 1); imagesc(map_phase); axis image; colormap(MapAnsColorMap);
    hold on; plot(plot_x, plot_y, '.k'); colorbar; title('Phase map'); hold off;
    subplot(2, 1, 2); imagesc(map_strength); axis image; 
    hold on; plot(plot_x, plot_y, '.k'); colorbar; title('AMP map'); hold off;
   % figure; 
    
case 'line'
    [map_strength, map_phase, file_name]= read_map('Please select the map', filter_size); %smoothing
    map_phase = map_phase/360 *deg_per_stimcycle; % deg in real distance
    hndl = figure; 
    imagesc(map_phase); axis image; colormap(MapAnsColorMap); colorbar; hold on; title('Phase map'); 
    disp('Please draw a line with mouse, 2 points');
    [rec_x, rec_y] = ginput(2);
    rec_x = round(rec_x); rec_y = round(rec_y);
    line_slope = diff(rec_y)/diff(rec_x);

    if abs(line_slope) <=1 
        if diff(rec_x)<0
            plot_x = [rec_x(1) : -1 : rec_x(2)]';
        else
            plot_x = [rec_x(1) : rec_x(2)]'; 
        end
        plot_y = round(rec_y(1) + (plot_x-rec_x(1))*line_slope);
   else
        if diff(rec_y)<0
            plot_y = [rec_y(1) : -1 : rec_y(2)]';
        else
            plot_y = [rec_y(1) : rec_y(2)]';
        end
        plot_x = round(rec_x(1) + (plot_y-rec_y(1))/line_slope);
    end
    plot(plot_x, plot_y, '.r'); hold off;

end

% To find all the points 
[ydim, xdim] = size(map_phase);
phase_plot = map_phase((plot_x-1)*ydim + plot_y); 
no_points = length(plot_x);

dist_from_first = sqrt((plot_x-plot_x(1)).*(plot_x-plot_x(1)) ...
    + (plot_y-plot_y(1)).*(plot_y-plot_y(1)))*mm_per_pixel; % in mm;
[dist_from_first, sort_index] = sort(dist_from_first);
phase_plot = phase_plot(sort_index);
plot_x = plot_x(sort_index);
plot_y = plot_y(sort_index);

hndl2 = figure(hndl+1);
plot(dist_from_first, phase_plot, '-k'); hold on;
% axis([-0.1 1.2 -10 50]);
xlabel('distance(mm)'); ylabel('phase(deg)');

dist_ends = dist_from_first(no_points) - dist_from_first(1);
phase_diff_ends = phase_plot(no_points) - phase_plot(1);
slope_ends = phase_diff_ends/dist_ends;

note_text = sprintf('%s:', file_name);
%print notes on the figure 
figure(hndl2);
axes('Position',[0 0 1 1],'Visible','off');
text(.025,0.95,note_text,'FontSize',10);

% convert to distance to the point of 0 phase
zero_phase = find(abs(phase_plot) == min(abs(phase_plot)));
zero_index = ceil(length(zero_phase)/2); % to choose the median point
dist_from_zero_phase = dist_from_first - dist_from_first(zero_phase(zero_index));

% % to calculate magnification factor between -20 deg and 20 deg on the map
m20_phase = find(abs(phase_plot+20) == min(abs(phase_plot+20)));
m20_index = ceil(length(m20_phase)/2);
p20_phase = find(abs(phase_plot-20) == min(abs(phase_plot-20)));
p20_index = ceil(length(p20_phase)/2);

phase_diff_20 = phase_plot(p20_phase(p20_index)) - phase_plot(m20_phase(m20_index));
dist_diff_20 = sqrt((plot_x(p20_phase(p20_index))-plot_x(m20_phase(m20_index)))^2 ...
    + (plot_y(p20_phase(p20_index))-plot_y(m20_phase(m20_index)))^2)*mm_per_pixel;
mag_factor = phase_diff_20 / dist_diff_20;

IMAPhase = struct ('MapName', file_name, ...
    'Distance', dist_from_zero_phase, 'Distance_first', dist_from_first,... 
    'Phase', phase_plot, 'MagFactor', mag_factor, ...
    'TemplateMethod', TemplateMethod, 'TemplateValue', TemplateValue, ...
    'Filter', filter_size, ...
    'Coordinates', [plot_x(m20_phase(m20_index)) plot_y(m20_phase(m20_index)) ...
        plot_x(p20_phase(p20_index)) plot_y(p20_phase(p20_index))]);

% coordinates: of the two points onthe contour, used to calculated maps
% tilt.


%% new codes for plotting response magnitude along the contour/line
%% 08/19/06
%% for this purpose, we need to extend the contour line, in order to plot
%% the smaller amplitude 
%%first, to linear fit the contour line
P = polyfit(plot_x, plot_y, 1);
x_index = [1:xdim];
y_index = round(P(1)*x_index +P(2));
out_bound = find(y_index<1 | y_index>ydim);
x_index(out_bound) = []; 
y_index(out_bound) = []; 
figure; set(gca, 'YDir', 'reverse'); hold on; axis([1 512 1 512]); 
plot(plot_x, plot_y, '.r'); 
plot(x_index, y_index, '-k'); 

amp_plot = map_strength((x_index-1)*ydim + y_index);
peak_index = find(amp_plot == max(amp_plot));
amp_plot = amp_plot ./amp_plot(peak_index(1)); %normalize to the peak
dist_from_first = sqrt((x_index-x_index(1)).^2 + (y_index-y_index(1)).^2)*mm_per_pixel;
figure; plot(dist_from_first, amp_plot, '-k'); hold on;
xlabel('distance(mm)'); ylabel('magnitude');

IMAamp = struct ('MapName', file_name, 'PlotX', x_index, 'PlotY', y_index, ...
    'Amp', amp_plot, 'Distance', dist_from_first);
%% 08/19/06------------------------------
