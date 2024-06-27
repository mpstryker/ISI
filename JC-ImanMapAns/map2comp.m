function map2comp(phase_or_amp)
% to compare the phase/amplitude of 2 maps
% Jianhua Cang, last modified 06-21-2006

global mm_per_pixel deg_per_stimcycle filter_size;
global TemplateMethod TemplateValue;
global MapAnsColorMap;

%% read in the existed template or make a new template

[response non_response template_file_name template_ydim template_xdim] = ...
        maketemplate(TemplateMethod, TemplateValue, filter_size);

%% the first map
    [map1_strength_filtered, map1_phase_filtered, file_name_1, map1_strength, map1_phase]= ...
    read_map('the first map :', filter_size);

%% the second map
[map2_strength_filtered, map2_phase_filtered, file_name_2, map2_strength, map2_phase]= ...
    read_map('the second map:', filter_size);

map1_phase = map1_phase/360 *deg_per_stimcycle; % deg in real distance
map2_phase = map2_phase/360 *deg_per_stimcycle; % deg in real distance
map1_phase_filtered = map1_phase_filtered/360 *deg_per_stimcycle; % deg in real distance
map2_phase_filtered = map2_phase_filtered/360 *deg_per_stimcycle; % deg in real distance

[ydim1, xdim1] = size(map1_strength);
[ydim2, xdim2] = size(map2_strength);
if(ydim1~=template_ydim)|(xdim1~=template_xdim)|(ydim2~=template_ydim)|(xdim2~=template_xdim)
    errordlg('map dimensions do not match!');
    return;
end

amp_ratio_map = (map2_strength - map1_strength)./(map2_strength + map1_strength);
amp_ratio_map_filtered = (map2_strength_filtered - map1_strength_filtered)...
        ./(map2_strength_filtered + map1_strength_filtered);
    
% amp_diff_map = map2_strength_filtered - map1_strength_filtered;
% amp_R_map = log10(map2_strength./map1_strength); % R distributions of individual pixels

switch phase_or_amp
    
case 'phase_comp'
    phase_diff_map = map2_phase - map1_phase; 
    phase_diff_map_filtered = map2_phase_filtered - map1_phase_filtered; 
    phase_diff_map(non_response) = 0;
    phase_diff_map_filtered(non_response) = 0;
    % average_shift = mean(phase_diff_map (response))
    hndl = figure; 
    subplot(2, 2, 1); 
    imagesc(map1_phase); colorbar; axis image; title(file_name_1);
    %imagesc(map1_phase_filtered); colorbar; axis image; title(file_name_1);
    subplot(2, 2, 2); 
    imagesc(map2_phase); colorbar; axis image; title(file_name_2);
    % imagesc(map2_phase_filtered); colorbar; axis image; title(file_name_2);
    subplot(2, 2, 3); 
    imagesc(phase_diff_map); axis image; colorbar; title('Phase Diff 2-1');
    subplot(2, 2, 4); 
    hist(phase_diff_map(response), [-40:40]); title('hist(Phase Diff)');
    set(gca,'XLim',[-40 40]);
    figure(hndl+1); imagesc(phase_diff_map_filtered); 
    axis image; colorbar; title('Phase Diff 2-1 Filtered');
    figure(hndl+2); hold on;
    plot(map2_phase_filtered(response), phase_diff_map_filtered(response), 'k.');
%     SumPhaseDiff = full(sparse(1, round(map2_phase_filtered(response))+10, phase_diff_map_filtered(response)));
%         % create a sparse matrix to calcualte sum, shift of 10 to take care
%         % of 0 and negative phase of map2
%     SumPhaseCount = full(sparse(1, round(map2_phase_filtered(response))+10, 1));
%         % index of how many points were included
%     PhaseDiffvsPhase = SumPhaseDiff./SumPhaseCount;
%     plot([-9:length(SumPhaseDiff)-10], PhaseDiffvsPhase, 'r*-');
        
case 'amp_comp'
    map1_strength(non_response) = 0;
    map2_strength(non_response) = 0;
    amp_ratio_map(non_response) = -1;
    amp_ratio_map_filtered(non_response) = -1;
    hndl = figure;      
    % plot the non-filtered maps
    color_range = [-max(max(map1_strength_filtered(:)), max(map2_strength_filtered(:))) ...
            -min(min(map1_strength_filtered(:)), min(map2_strength_filtered(:)))];
    % color_range = [-2.4 0];
    subplot(2,2,1); imagesc(-map1_strength, color_range); colorbar; axis image; title('map 1');
    subplot(2,2,2); imagesc(-map2_strength, color_range); colorbar; axis image; title('map 2');
    subplot(2,2,3); imagesc(-amp_ratio_map, [-0.5 0.5]); colorbar; axis image; title('nonfiltered (2-1)/(2+1)');
            % for making figures, grayer-stronger
    subplot(2,2,4); % histogram of the unfiltered OD map
    hist(amp_ratio_map(response), [-1:0.02:1]); title('nonfiltered (2-1)/(2+1)'); 
    set(gca,'XLim',[-1 1]); colormap('gray');
      
    figure(hndl+1); 
    % imagesc(amp_ratio_map, [-0.5159 0.5]); colorbar; axis image; title('non filtered (2-1)/(2+1)');
    imagesc(amp_ratio_map_filtered, [-0.5159 0.5]); colorbar; axis image; title('filtered (2-1)/(2+1)');
    newcolor= [[1 1 1]; colormap];
    colormap(newcolor);
    
    % measurements of the ratio map
    response_pixels = length(response);
    ratio_map_area = response_pixels * mm_per_pixel * mm_per_pixel; % in mm2
    ratio_map_volume = sum(amp_ratio_map(response))* mm_per_pixel * mm_per_pixel;
    ratio_map_mean_height = sum(amp_ratio_map(response))/response_pixels; %unfiltered
    
    %print notes in command window
    disp(sprintf('Max response: %s %3.2f\n%s %3.2f and ratio(I/C): %3.2f', ...
        file_name_1, max(map1_strength_filtered(response)), file_name_2, ...
        max(map2_strength_filtered(response)), ...
        max(map1_strength_filtered(response))/max(map2_strength_filtered(response))));
    %disp(sprintf('Minimum of ratio map is %3.2f', amp_min_ratio));
    disp(sprintf('Area of ratio map is %3.2f', ratio_map_area));
    disp(sprintf('Volume of ratio map is %3.2f', ratio_map_volume));
    disp(sprintf('Average Height of ratio map is %3.2f', ratio_map_mean_height));
   
    %print notes on the figure 
    note_text = sprintf('%s and %s, C-I/C+I: %3.2f \n Template: %s', ...
        file_name_1, file_name_2, ratio_map_mean_height, template_file_name);
    figure(hndl);
    axes('Position',[0 0 1 1],'Visible','off');
    text(.025,0.035,note_text,'FontSize',10);
    
    figure(hndl+1);
    axes('Position',[0 0 1 1],'Visible','off');
    text(.025,0.025,note_text,'FontSize',10);
end

