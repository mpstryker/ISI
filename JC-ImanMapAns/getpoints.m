function [xs,ys, linehandle] = getpoints(axishandle,dospline)

%============================================================================
% [xs,ys, linehandle] = getpoints(axishandle,dospline)
% select an area interactively with mouse clicking a polygon
% right click to close the selection
% Jianhua Cang, Dec-12-2005, with original code from Matt Caywood
%============================================================================

%============================================================================
% Find parent figure for the argument axishandle
%============================================================================
axes(axishandle);
figure(get(axishandle, 'Parent'));
%===========================================================================
% Change pointer shape
%===========================================================================
% oldpointershape = get(gcf,'Pointer');
% 
% ptrc =  ones(16)+1;
% ptrc( 1, :) = 1;
% ptrc(16, :) = 1;
% ptrc(: , 1) = 1;
% ptrc(: ,16) = 1;
% ptrc(1:4,8:9) = 1;
% ptrc(8:9,1:4) = 1;
% ptrc(13:16, 8:9 ) = 1;
% ptrc( 8:9 ,13:16) = 1;
% ptrc(5:12,5:12) = NaN;
% set(gcf,'Pointer',              'custom',...
%        'PointerShapeCData',    ptrc,...
%         'PointerShapeHotSpot', [8 8]);

%===========================================================================
% Prepare for interactive collection of ROI boundary points
%===========================================================================
hold on
pointhandles = [];
xpts = [];
ypts = [];
outlinehandle= [];
n = 0;
but = 1;
BUTN = 0;
KEYB = 1;
done =0;

%===========================================================================
% Loop until right hand mouse button or keayboard is pressed
%===========================================================================
while ~done;
 %===========================================================================
 % Analyze each buttonpressed event
 %===========================================================================
 keyb_or_butn = waitforbuttonpress;
 if keyb_or_butn == BUTN;
   currpt = get(axishandle, 'CurrentPoint');
   seltype = get(gcf,'SelectionType');
   switch seltype
   case 'normal',
     but = 1;
   case 'alt',
     but = 2;
   otherwise,
     but = 2;
   end;
 elseif keyb_or_butn == KEYB
   but = 2;
 end;
 %===========================================================================
 % Get coordinates of the last buttonpressed event
 %===========================================================================
 xi = currpt(2,1);
 yi = currpt(2,2);
 %===========================================================================
 % Start a spline throught the points or
 % update the line through the points with a new spline
 %===========================================================================
 if but ==1
   if ~isempty(outlinehandle)
      delete(outlinehandle);
   end;
   pointhandles(n+1) = plot(xi,yi,'wo','MarkerSize',4,'LineWidth',2);
       n = n+1;
       xpts(n,1) = xi;
       ypts(n,1) = yi;
       %===========================================================================
       % Draw a spline line through the points
   %===========================================================================
       if n > 1
         t = 1:n;
         ts = 1: 0.1 : n;
         if (dospline)
         xs = spline(t, xpts, ts);
         ys = spline(t, ypts, ts);
         outlinehandle = plot(xs,ys,'r-');
     else
         outlinehandle = line(xpts,ypts,'Color',[1 1 1],'LineWidth',2);
     end
       end;
 elseif but > 1
     %===========================================================================
         % Exit for right hand mouse button or keyboard input
     %===========================================================================
     done = 1;
 end;
end;

%===========================================================================
% Add first point to the end of the vector for spline
%===========================================================================
xpts(n+1,1) = xpts(1,1);
ypts(n+1,1) = ypts(1,1);

%===========================================================================
% (re)draw the final spline
%===========================================================================
if ~ isempty(outlinehandle)
   delete(outlinehandle);
end;

if (dospline)
   t = 1:n+1;
   ts = 1: 0.25 : n+1;
   xs = spline(t, xpts, ts);
   ys = spline(t, ypts, ts);

   linehandle = plot(xs,ys,'r-');
else
   linehandle = line(xpts,ypts,'Color',[1 1 1],'LineWidth',2);
   xs = xpts;
   ys = ypts;
end

drawnow;
%===========================================================================
% Delete the point markers
%===========================================================================
if ~isempty(pointhandles)
   delete(pointhandles)
end;
%===========================================================================
% Reset pointershape
%===========================================================================
% set(gcf,'Pointer',oldpointershape);
