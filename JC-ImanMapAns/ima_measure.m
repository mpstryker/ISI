function IMA_measure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab code for calculating response area
% Two methods to set threshold: X fraction of peak response 
%                            or X folds of background
% Jianhua Cang, 04-16-04
% last updated 06-21-06
%%%%%%%%%%%%%%%%%%%%%%%%%%

global deg_per_stimcycle mm_per_pixel filter_size;
global TemplateMethod TemplateValue;

%% display these value to be sure the program read into correct parameters
deg_per_stimcycle 
mm_per_pixel 
filter_size
TemplateMethod 
TemplateValue

[Amap_strength, Amap_phase, file_name] = read_map('Please select the map', filter_size);
Amap_phase = Amap_phase/360 *deg_per_stimcycle; 
[ydim, xdim] = size(Amap_strength);
hndl = figure;  subplot(2,2,1);
imagesc(Amap_strength); axis image; colorbar; 
[X, Y] = meshgrid(1:xdim, 1:ydim);
axishandle = gca;

if (TemplateMethod == 3)
    bkgrd_select = questdlg('Please select a background','bkgrd_select','Yes','OK','OK');
    disp('select a polygon with mouse, right click to close');
    [xPoly, yPoly, linehandle] = getpoints(axishandle,0);
    pts_inside = inpolygon(X, Y, xPoly, yPoly);
    bkgrd = mean(Amap_strength(pts_inside));
    threshold_height = TemplateValue * bkgrd;
end

region_select = questdlg('Select a region for analysis?',...
    'Region_select','Yes','No','No');

if strcmp(region_select,'Yes')
     disp('select a polygon with mouse, right click to close');
     [xPoly, yPoly, linehandle] = getpoints(axishandle,0);
     h = waitbar(0,'Analyzing ...');
     pts_inside = inpolygon(X, Y, xPoly, yPoly);            
     h = waitbar(0.8,h);
     Amap_strength(~pts_inside) = 0;  
     waitbar(1, h); close(h);
     figure(hndl); subplot(2,2,3); 
     imagesc(Amap_strength); axis image; colorbar;title('After select region'); 
end
max_height = max(Amap_strength(:));

if (TemplateMethod == 2)
    threshold_height=TemplateValue * max_height;
end
    
non_response = find(Amap_strength < threshold_height);
Amap_strength(non_response) = 0;
response = find(Amap_strength >= threshold_height);
area = length(response) * mm_per_pixel * mm_per_pixel; % in mm2
% volume = sum(filtered_map(response)-threshold_height)* mm_per_pixel * mm_per_pixel

figure(hndl);
subplot(2, 2, 2); ;
imagesc(Amap_strength); axis image; colorbar;title('After thresholding')
subplot(2, 2, 4); 
imagesc(Amap_phase); axis image; colorbar;title('Phase map');

note_text = sprintf('%s: Max: %3.2f; Area: %3.2f mm2', ...
        file_name, max_height, area);
%print notes on the figure 
figure(hndl);
axes('Position',[0 0 1 1],'Visible','off');
text(.025,0.95,note_text,'FontSize',10);