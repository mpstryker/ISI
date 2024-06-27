function OverLapPercent = OverLap2Circles(R1, R2, D)

%% to calculate the overlap between two circles
%% the radii are symbolized as R1 and R2, and the distance between the centers is D 

if (R1+R2)<=D %% no overlap
    OverLapPercent = 0;
    return
end

if (R1+D)<=R2 %% R1 is inside R2
    OverLapPercent = (R1*R1)/(R2*R2);
    return
end

if (R2+D)<=R1 %% R2 is inside R1
    OverLapPercent = (R2*R2)/(R1*R1);
    return
end

% we know the three sides of the triangle between the circle
% centers and one intersection. The angle between the the line
% connecting the centers and a line drawn from the center of first
% circle to one intersection is Phi1, and is given by this formula:
 
 Phi1 = acos((D^2 + R1^2 - R2^2)./(2 * R1 * D));
 
%  If we consider the area of the overlap to be considered
% as two adjacent circular "segments", then the area of the
% segment farthest from the center of the first circle is
% given by these formulas:

Theta1 = 2*Phi1;
Area1 = ( R1^2 * (Theta1 - sin(Theta1)) ) / 2;

% Do this for the other segment

Phi2 = acos((D^2 + R2^2 - R1^2)./(2 * R2 * D));
 
%  If we consider the area of the overlap to be considered
% as two adjacent circular "segments", then the area of the
% segment farthest from the center of the first circle is
% given by these formulas:

Theta2 = 2*Phi2;
Area2 = ( R2^2 * (Theta2 - sin(Theta2)) ) / 2;

OverlapArea = Area1+Area2;

OverLapPercent = OverlapArea / (pi*R1*R1 + pi*R2*R2 -OverlapArea );