% Matlab code for reading maps genereted by 
% iman -s -D -N32 -T1 (which can be viewed binan -i filename
% Jianhua Cang, 10-18-04

function varargout = read_movie(varargin)

global xdim ydim;
global Nbin;

% PromptText = 'Select Data File';
PromptText = varargin{1};;

% to filter or not?
% if nargin == 1
%     filter_or_not = false;
% elseif nargin == 2
%     filter_or_not = true; 
%     filter_size = varargin{2};
%     filter_matrix = ones(filter_size)/(filter_size*filter_size);
% end 

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
Nbin = fread(fid,1,'int32');
array_size = xdim*ydim;

h = waitbar(0,'reading maps ...');
for i = 1: Nbin
    temp = fread(fid,array_size,'double');
    strength_matrix(:, :, i) = reshape(temp, xdim, ydim)'; 
    waitbar(i/Nbin,h)
end
close(h);
fclose(fid);

% hndl = figure;
% for i = 1: Nbin
%     colormap('gray');
%     subplot(4, ceil(Nbin/4), i);
%     imagesc(strength_matrix(:, :, i)); title(sprintf('%d', i));
% end

% note_text = sprintf('%s:', data_filename);
% %print notes on the figure 
% figure (hndl); 
% axes('Position',[0 0 1 1],'Visible','off');
% text(.025,0.95,note_text,'FontSize',10);


disp(sprintf('xdim is %d', xdim));
disp(sprintf('ydim is %d', ydim));
disp(sprintf('number of frames is %d', Nbin));


varargout(1) = {strength_matrix};
varargout(2) = {Nbin};
varargout(3) = {data_filename};