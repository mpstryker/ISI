function img3 = NearScreen(img, cp_azdeg, cp_eldeg, cp_distcm, cpx_cm, cpy_cm, degPerPix, sizeX, sizeY, pixelsPerCm, sizeXcm, sizeYcm)

[cp_azdeg cp_eldeg cp_distcm cpx_cm cpy_cm degPerPix sizeX sizeY pixelsPerCm sizeXcm sizeYcm]

% cp_azdeg=0; cp_eldeg=0; cp_distcm=12.5; 
% cpx_cm=20; cpy_cm=10; pixelsPerCm=20;
% sizeXcm=40;
% sizeYcm=30;

% orient=90;
% freq=.5; 
% speed=30;
% contrast=.95;
% length=0;
% position=0;
% duration=.1;
% sizeX=800;
% sizeY=600; 
% frameRate=60; 
% black=0; 
% white=255;
% sizeLut=256;
% degPerPix=sizeXcm*atan(1/cp_distcm)*(180/pi)/sizeX;
% disp 'img '; tic
% [img lut]= generateBars_lut(orient,freq,speed,contrast,length, position, duration, degPerPix,sizeX,sizeY, frameRate, black, white,sizeLut);
% im8 = uint8(img);
% image(im8, 'CDataMapping', 'scaled');
% figure;
% image(img, 'CDataMapping', 'scaled');
% toc
disp 'img2 '; tic
%img2=257*ones(600,800);
img2=zeros(sizeY,sizeX);
for x=[1:.25:sizeX] 
    for y=[1:1:sizeY]
       az_deg =   cp_azdeg + ((((x/sizeX)*sizeXcm) - cpx_cm)*(pixelsPerCm*degPerPix));
       elev_deg = cp_eldeg + ((((y/sizeY)*sizeYcm) - cpy_cm)*(pixelsPerCm*degPerPix));
        
        [Pix_x Pix_y] = pt2screen(az_deg, elev_deg, cp_azdeg, cp_eldeg, cp_distcm, cpx_cm, cpy_cm, pixelsPerCm);
        if ((Pix_y > 0) & (Pix_y <= sizeY) & (Pix_x > 0) & (Pix_x <= sizeX))
            img2(Pix_y,Pix_x) = img(round(y),round(x));  %?? sizeY-Pix_y ??
        end
    end
end
figure;
image(img2, 'CDataMapping', 'scaled');
toc
disp 'filter '; tic
img3=img2;
ind = find(img3 == 0);
for k=[5:size(ind,1)-5]
    kk=ind(k);
    img3(kk)=max(img2(kk-3:kk+3));
end

%img3 = medfilt2(img2,[5 5]);
figure;
image(img3, 'CDataMapping', 'scaled');
toc
% 

