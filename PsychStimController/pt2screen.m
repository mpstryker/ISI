function [Pix_x Pix_y] = pt2screen(az_deg, elev_deg, ...
    cp_azdeg, cp_eldeg, cp_distcm, cpx_cm, cpy_cm, pixelsPerCm)
% function to convert a point specified as absolute azimut and elevation to
% pixels on a screen.  
% Screen position is specified wrt its closest point by the the 
% azimuth and elevation and distance of that point in deg and cm, 
% plus the location of that closest point on the screen in cm.
% pixelsPerCm = 10000/pixelwidth_um = screenWidth_pixels / screenWidth_cm;

% cpaed is closest point azimuth elevation and distance wrt straight ahead
acp = deg2rad(cp_azdeg);
ecp = deg2rad(cp_eldeg);
dcp = cp_distcm;

% cpxyz is closest point on sphere in cartesian coordinates, 
% where x is to right (East), y is up (North), and z is straight ahead
cpx = dcp*(cos(ecp) * sin (acp));
cpy  = dcp*(sin(ecp));
cpz = dcp*(cos(ecp) * cos (acp));
%[cpx cpy cpz]
% <a,e> is target location 
a = deg2rad(az_deg); 
e = deg2rad(elev_deg);
% Screen Plane is: cpx*x + cpy*y + cpz*z - dcp = 0;
%P1 = <0,0,0> = <x1,y1,z1>
%P2 = <(cos(e) * sin(a)), sin(e), (cos(e) * cos(a))> = <x2,y2,z2>
% intersection P = P1 + u(P2-P1)
% Since P1 = <0,0,0>, then intersection = u*P2
u = (dcp)/( (cpx*(cos(e)*sin(a))) + (cpy*sin(e)) + (cpz*(cos(e)*cos(a))) );
% Intersection point of target line with screen sc: <scx,scy,scz> in cartesian coordinates
scx = dcp*u*(cos(e)*sin(a));
scy = dcp*u*(sin(e));

scz = dcp*u*(cos(e)*cos(a));
%[scx scy scz]
% Intersection point s <sx,sy,sz> in Screen coordinates wrt cp is projection
% onto screen plane
scrx = (scx - cpx)/cos(acp);
scry = (scy - cpy)/cos(ecp);

Pix_x = fix((cpx_cm + scrx)* pixelsPerCm);
Pix_y = fix((cpy_cm + scry)* pixelsPerCm);
end
