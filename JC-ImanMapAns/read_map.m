% Matlab code for reading maps genereted by mapans
% read in map and filter
% Jianhua Cang, 09-24-04

function varargout = read_map(varargin)

%[strength_matrix_filtered, phase_matrix_filtered, 
% data_filename, strength_matrix, phase_matrix]=read_map(varargin)

global xdim ydim;

% PromptText = 'Select Data File';
PromptText = varargin{1};;

% to filter or not?
if nargin == 1
    filter_or_not = false;
elseif nargin == 2
    filter_or_not = true; 
    filter_size = varargin{2};
    filter_matrix = ones(filter_size)/(filter_size*filter_size);
end 

[data_filename, data_pathname] = uigetfile('*.*', PromptText);
if isequal(data_filename,0) | isequal(data_pathname,0)
   disp('File not found')
   return;
else
   disp(['You selected ', fullfile(data_pathname, data_filename)])
end

cd(data_pathname);
% data_filename = [data_pathname data_filename]
fid = fopen(data_filename, 'r');
type = fread(fid,1,'int32');
xdim = fread(fid,1,'int16');
ydim = fread(fid,1,'int16');
array_size = xdim*ydim;
x_array = fread(fid,array_size,'double');
y_array = fread(fid,array_size,'double');
fclose(fid);

disp(sprintf('xdim is %d', xdim));
disp(sprintf('ydim is %d', ydim));

% change array to matrix
x_matrix = reshape(x_array, xdim, ydim)';  
y_matrix = reshape(y_array, xdim, ydim)';
xy_matrix = x_matrix + i*y_matrix;

phase_matrix = angle(xy_matrix);
phase_matrix = round(phase_matrix*180/pi);
strength_matrix = abs(xy_matrix);

phase_matrix_filtered = phase_matrix;
strength_matrix_filtered = strength_matrix;

if (filter_or_not)
    xy_matrix_filtered = filter2(filter_matrix, xy_matrix);
    phase_matrix_filtered = angle(xy_matrix_filtered);
    phase_matrix_filtered = round(phase_matrix_filtered*180/pi);
    strength_matrix_filtered = abs(xy_matrix_filtered);
end

varargout(1) = {strength_matrix_filtered};
varargout(2) = {phase_matrix_filtered};
varargout(3) = {data_filename};
varargout(4) = {strength_matrix};
varargout(5) = {phase_matrix};

% warning off MATLAB:divideByZero;
