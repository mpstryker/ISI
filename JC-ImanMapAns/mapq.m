
function varargout = mapQ(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab code for quantifying map quality
% read in raw map without filtering 
% then filter, calculate difference between filtered map and raw map
% Jianhua Cang, 04-16-04
% last updated 09-30-05
%%%%%%%%%%%%%%%%%%%%%%%%%%
global mm_per_pixel deg_per_stimcycle;
global TemplateMethod TemplateValue;
global MapAnsColorMap;

HistRange = [-90:90];

filter_size = 5; 
    % default filter size, with odd number, the filtered matrix does not
    % shift position, checked on 11-25-04
if nargin == 1
    filter_size = varargin{1};
end 

[map_strength_filtered, map_phase_filtered, file_name, map_strength, map_phase]= ...
    read_map('Please select the target map:' , filter_size);
map_phase = map_phase/360 *deg_per_stimcycle; % deg in real distance
map_phase_filtered = map_phase_filtered/360 *deg_per_stimcycle; % deg in real distance
[ydim, xdim] = size(map_phase);

if (TemplateMethod ~= 6) % NOT local analysis
    [ROI non_ROI template_file_name template_ydim template_xdim] = maketemplate(TemplateMethod, TemplateValue, filter_size);
    
    diff_from_filtered = map_phase - map_phase_filtered;
    hdnl1 = figure; %colormap(gray);
    subplot(2,2,1); imagesc(map_phase); colormap(MapAnsColorMap); colorbar; title('raw map');
    subplot(2,2,3); imagesc(diff_from_filtered); colormap(MapAnsColorMap);  
    colorbar; title('Phase difference: raw-filtered')
%     diff_from_filtered(non_ROI) = 0; 
%     subplot(2,2,2); imagesc(diff_from_filtered); colorbar; title('diff filtered')
    subplot(2,2,4); 
    bar(HistRange, hist(diff_from_filtered(ROI), HistRange), 'k'); 
    title('ROI, diff filtered');
    subplot(2,2,2); 
    bar(HistRange, hist(diff_from_filtered(non_ROI), HistRange), 'k'); 
    title('outside ROI, diff filtered');
      
    % stats
    % mean(abs(diff_phase(ROI)))
    std_diff_filtered =std(diff_from_filtered(ROI));
    note_text_1 = sprintf('%s: Filter Size, %d', file_name, filter_size);
    note_text_2 = sprintf('ROI: std(rawmap-filtered): %3.2f', ...
        std_diff_filtered);
    
%     %%%%%%%07-27-05%%%%%%%%
%     % plot difference map
%     figure; color('gray');
%     imagesc(abs(diff_from_filtered), [0 10]); colorbar; title('Phase difference: abs(raw-filtered)')
    
else % local analysis 
    %% needs work to move this into maketemplate
    %%%%%%% ask the user input position of interest and area surrouding
    if (TemplateValue <=1) 
        errordlg('the diameter is too small');
        return;
    end
    ROI_diameter = TemplateValue; % in pixels
    hndl2 = figure; imagesc(map_phase); colormap(MapAnsColorMap); colorbar; 
    hold on; title('raw map');
    prompt = {'Enter x:','Enter y:'};
    dlg_title = 'Input for ROI coordinates';
    temp = inputdlg(prompt,dlg_title);
    temp =cell2mat(temp);
    ROI_x = str2num(temp(1, :));    ROI_y = str2num(temp(2, :));
    figure(hndl2);
    plot(ROI_x, ROI_y, 'xk'); set(gca,'YDir','reverse');
    
    % creat position index for calculating distance: 
    % colum_index and row_index
    x_index = [1:xdim]; 
    y_index = [1:ydim]';
    ro_index = y_index*ones(1, xdim); %y index
    co_index = ones(ydim, 1)*x_index; % x index
    distance_matrix = sqrt((ro_index-ROI_y).*(ro_index-ROI_y) ...
       +(co_index-ROI_x).*(co_index-ROI_x));
    ROI = find(distance_matrix <= ROI_diameter/2);
    non_ROI = find(distance_matrix > ROI_diameter/2);
    %save the template in a struture
    save_template = questdlg('Do you want to save the template?',...
        'Save Template','Yes','No','Yes');
    if strcmp(save_template,'Yes')
        Save_File_template = [file_name(1:6) 'template.mat'];
        TemplateMap = struct ('Name', file_name, 'ydim', ydim, 'xdim', xdim', ...
        'response', ROI, 'non_response', non_ROI, ...
        'center_x', ROI_x, 'center_y', ROI_y);
        [filename, pathname] = uiputfile( '*.mat', 'Save Workspace as', Save_File_template);
        save(fullfile(pathname, filename), 'TemplateMap');  
        
    end

    diff_from_center = map_phase-map_phase(ROI_y, ROI_x);
    diff_from_filtered = map_phase - map_phase_filtered;
    
    hdnl1 = figure; %colormap(gray);
    subplot(2,3,1); imagesc(map_phase); 
    colormap(MapAnsColorMap); colorbar; hold on; title('raw map');
    plot(ROI_x, ROI_y, 'xk'); set(gca,'YDir','reverse');
    subplot(2,3,4); imagesc(diff_from_filtered); 
    colormap(MapAnsColorMap); colorbar; title('Phase difference: raw-filtered')
    diff_from_filtered(non_ROI) = 0; 
    subplot(2,3,2); imagesc(diff_from_filtered); 
    colormap(MapAnsColorMap); colorbar; title('diff filtered')
    subplot(2,3,5); hist(diff_from_filtered(ROI), [-20:20]); title('ROI, diff filtered');
    diff_from_center(non_ROI) = 0; 
    subplot(2,3,3); imagesc(diff_from_center); 
    colormap(MapAnsColorMap); colorbar; title('diff center')
    subplot(2,3,6); hist(diff_from_center(ROI), [-20:20]); title('ROI, diff center');
     
    % stats
    % mean(abs(diff_phase(ROI)))
    std_diff_center =std(diff_from_center(ROI));
    std_diff_filtered =std(diff_from_filtered(ROI));
    note_text_1 = sprintf('%s: Filter Size, %d', file_name, filter_size);
    note_text_2 = sprintf('ROI: std(rawmap-filtered): %3.2f; std(rawmap-center): %3.2f', ...
        std_diff_filtered, std_diff_center);
end
   
%print notes on the figure 

figure(hdnl1); 
axes('Position',[0 0 1 1],'Visible','off');
text(.025,0.05,note_text_1,'FontSize',10);
axes('Position',[0 0 1 1],'Visible','off');
text(.025,0.025,note_text_2,'FontSize',10);


