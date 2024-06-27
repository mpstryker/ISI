function [response, non_response, file_name, ydim, xdim] = maketemplate(Method, Value, filter_size)

% Make a template using one of the following method
% Method = 1: exsited template
% Method = 2: percentage of peak response as threshold, set by Value
% Method = 3: fold of backgroud response as threshold, set by Value
% Method = 4: Number of points with strongest response, set by Value
% Method = 5: TIF file from IDL
% Method = 6: local coordinates
% last modified, 12-12-05

global MapAnsColorMap;
global deg_per_stimcycle;

switch Method
    case 1 %% exsited template
        [data_filename, data_pathname] = uigetfile('*.mat', 'select the template');
        load(fullfile(data_pathname, data_filename));
        file_name = TemplateMap.Name;
        ydim = TemplateMap.ydim;
        xdim = TemplateMap.xdim;
        response = TemplateMap.response;
        non_response = TemplateMap.non_response;
   
    case {2, 3, 4} 
        if (Method == 2 & Value >= 1)| ((Method == 3 & Value < 1))
            errordlg('the threshold is not appropriate for this method');
            return;
        end

        [map_strength_filtered, map_phase_filtered, file_name] = read_map('Select map for making template', filter_size);
        temp_map_strength = map_strength_filtered; % for later anlysis
        map_phase_filtered = map_phase_filtered/360 *deg_per_stimcycle; % deg in real distance
        [ydim, xdim] = size(map_strength_filtered);
        [X, Y] = meshgrid(1:xdim, 1:ydim);
        
        OutputFigure = figure;
        subplot(2, 2, 1); imagesc(map_strength_filtered);
        axis image; colormap('default'); colorbar; title('Original map for making template');
        subplot(2, 2, 3); imagesc(map_phase_filtered); 
        axis image; colorbar; title('Phase map');
        %%colormap(MapAnsColorMap); %% need quite some work to use 2 color maps
        %%on subplots of the same figure, maybe later
        
        InteractiveFigure = figure; imagesc(map_strength_filtered); 
        axis image; colormap('default'); colorbar; title('Original map for making template');
        axishandle = gca; 
           
        region_select = questdlg('Do you want to select a region?',...
                                 'Region_select','Yes','No','No');

        if strcmp(region_select,'Yes')
            figure(InteractiveFigure);
            title('select a polygon with mouse, right click to close');
            [xPoly, yPoly, linehandle] = getpoints(axishandle,0);
            h = waitbar(0,'Analyzing ...');
            pts_inside = inpolygon(X, Y, xPoly, yPoly);            
            h = waitbar(0.8,h);
            temp_map_strength(~pts_inside) = 0;  
            waitbar(1, h); close(h);
            figure(OutputFigure);
            subplot(2, 2, 2); imagesc(temp_map_strength); 
            axis image; colorbar; title('After select region');
        end

        if (Method == 2)
            max_height = max(map_strength_filtered(:));
            threshold_height=Value * max_height;
            Save_File_template = [file_name(1:6) 'templatePEAK.mat'];
        end

        if (Method == 3)
            figure(InteractiveFigure);
            bkgrd_select = questdlg('Please select a background?','bkgrd_select','Yes','OK','OK');            
%             disp('select a polygon with mouse, right click to close');
            [xPoly, yPoly, linehandle] = getpoints(axishandle,0);
            pts_inside = inpolygon(X, Y, xPoly, yPoly);
            bkgrd = mean(map_strength_filtered(pts_inside));
            threshold_height=Value * bkgrd;
            Save_File_template = [file_name(1:6) 'templateBKGRD.mat'];
          end
        
        if (Method == 4)
            temp = sort(map_strength_filtered(:));
            threshold_height=temp(numel(map_strength_filtered)-Value);
            Save_File_template = [file_name(1:6) 'templateSNP.mat'];
        end
        
        % to filter background    
        
        response = find(temp_map_strength >= threshold_height);
        non_response = find(temp_map_strength < threshold_height);
        temp_map_strength(non_response) = 0;
        figure(OutputFigure);
        subplot(2, 2, 4);
        imagesc(temp_map_strength); colormap('default');
        axis image; colorbar;title('The template');
       
               
    case 5       
        [data_filename, data_pathname] = uigetfile('*.tif', 'Select the TIF file');
        if isequal(data_filename,0) | isequal(data_pathname,0)
            errordlg('File not found'); return;
        end
        disp(['You selected ', fullfile(data_pathname, data_filename)]);
                cd(data_pathname);
        % data_filename = [data_pathname data_filename]
        idl_map = imread(data_filename);
        [ydim, xdim] = size(idl_map);
        hndl = figure; imagesc(idl_map); axis image; colormap('gray');
        response = find(idl_map);
        non_response = find(idl_map == 0);
        Save_File_template = [data_filename(1:6) 'IDLtemplate.mat'];
        
    case 6  %% local coordinates, 
        %%%%%%%%%%%%Note on 09-30-05, needs work
        %%%%%%% ask the user input position of interest and area surrouding
        %%%%%%% for now, use 50 pixels
        ROI_diameter = 50; % in pixels
        [map_strength_filtered, map_phase_filtered, file_name, map_strength, map_phase] ...
            = read_map('Select map for making template', filter_size);
        [ydim, xdim] = size(map_strength_filtered);
        hndl2 = figure; imagesc(map_phase); colormap(MapAnsColorMap); colorbar; 
        hold on; title('raw map');
        prompt = {'Enter x:','Enter y:'};
        dlg_title = 'Input for ROI coordinates';
        temp = inputdlg(prompt,dlg_title);
        temp =cell2mat(temp);
        ROI_x = str2num(temp(1, :));    ROI_y = str2num(temp(2, :));
        Value = temp;
        figure(hndl2); plot(ROI_x, ROI_y, 'xk'); set(gca,'YDir','reverse');
    
        % creat position index for calculating distance: 
        % colum_index and row_index
        x_index = [1:xdim]; 
        y_index = [1:ydim]';
        ro_index = y_index*ones(1, xdim); %y index
        co_index = ones(ydim, 1)*x_index; % x index
        distance_matrix = sqrt((ro_index-ROI_y).*(ro_index-ROI_y) ...
           +(co_index-ROI_x).*(co_index-ROI_x));
        response = find(distance_matrix <= ROI_diameter/2);
        non_response = find(distance_matrix > ROI_diameter/2);
        Save_File_template = [file_name(1:6) 'templateLocal.mat'];
end
        
%save the template in a struture
if (Method ~= 1)
    save_template = questdlg('Do you want to save the template?',...
        'Save Template','Yes','No','Yes');
    if strcmp(save_template,'Yes')
        TemplateMap = struct ('Name', file_name, 'ydim', ydim, 'xdim', xdim', ...
            'response', response, 'non_response', non_response, 'Method', Method, 'Value', Value);
        [filename, pathname] = uiputfile( '*.mat', 'Save Workspace as', Save_File_template);
        save(fullfile(pathname, filename), 'TemplateMap');  
    end
end