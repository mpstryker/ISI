function Contours = makecontour (varargin)

% Codes for plotting iso-phase contour
% For comparing map properties
% Last modified 09-30-2005
% Jianhua Cang
%%%%%%%%%%%%%%% input format:
%%%%%%%%%%%%%%%makecontour(isovalue1, isovalue2, isovalue3, ...)

global mm_per_pixel deg_per_stimcycle filter_size;
global TemplateMethod TemplateValue;
global MapAnsColorMap;

clear iso_value;
for i = 1:nargin
    iso_value(i) = varargin{i};
end

[map_strength, map_phase, file_name]= read_map('The map for drawing contour', filter_size);  % smoothed map
map_phase = map_phase/360 *deg_per_stimcycle; % deg in real distance
[ydim, xdim] = size(map_strength);

[response non_response template_file_name template_ydim template_xdim] = ...
        maketemplate(TemplateMethod, TemplateValue, filter_size);

if(ydim~=template_ydim)|(xdim~=template_xdim)
    errordlg('map dimensions do not match!');
    return;
end

map_strength(non_response) = 0;
figure; imagesc(map_phase); colormap(MapAnsColorMap); axis image; hold on; colorbar;

number_contours = length(iso_value);
clear iso_x iso_y;

% JC's way to plot contour before 10-12-04
% first to get the resolution of contour line, default 0 
% contour_res = abs(str2num(get(findobj('tag','contourresolutionedit'), 'string')));
% Contours = struct ('Name', file_name, 'Resolution', contour_res);
% for i = 1: number_contours 
%     [temp_y, temp_x] = find(map_phase >= (iso_value(i)-contour_res) & map_phase <= (iso_value(i)+contour_res) & map_strength);
%     plot(temp_x, temp_y, '.k');
%     Contours(i).isoValue = iso_value;
%     Contours(i).isoY = temp_y;
%     Contours(i).isoX = temp_x;
% end

% using matlab function contour
Contours = struct ('Name', file_name);
for i = 1: number_contours 
    [C, h]= contour(map_phase, [iso_value(i) iso_value(i)]);
    index = find(C(1, :) ~= iso_value(i));
    C= C(:, index); % remove 
    temp_x = round(C(1, :)); temp_y = round(C(2, :));
    index = find(temp_x); % remove all zeros
    temp_x = temp_x(index); temp_y = temp_y(index);
    index = find(temp_y); % remove all zeros
    temp_x = temp_x(index); temp_y = temp_y(index);
    
    index =find(map_strength(temp_y + (temp_x-1)*ydim)); % find within the template
    temp_x = temp_x(index); temp_y = temp_y(index);
    plot(temp_x, temp_y, '.k');

    region_select = questdlg('Do you want to select a region?',...
    'Region_select','Yes','No','No');
    if strcmp(region_select,'Yes')
        disp('select region with mouse, 2 points');
        [rec_x, rec_y] = ginput(2);
        rec_x = sort(round(rec_x));   rec_y = sort(round(rec_y));
        index_seleted = find(temp_x>rec_x(1) & temp_x<rec_x(2) & temp_y>rec_y(1) & temp_y<rec_y(2));
        temp_x = temp_x(index_seleted);  temp_y = temp_y(index_seleted);
        plot(temp_x, temp_y, '.-r');  
    end

    Contours(i).isoValue = iso_value;
    Contours(i).isoY = temp_y;
    Contours(i).isoX = temp_x;
    
end

save_dlg = questdlg('Do you want to save the contour coordinates?',...
        'Save','Yes','No','Yes');
if strcmp(save_dlg,'Yes')
    Save_File_Contour = [file_name 'Contour.mat'];
    [filename, pathname] = uiputfile( '*.mat', 'Save Workspace as', Save_File_Contour);
    save(fullfile(pathname, filename), 'Contours'); 
end
