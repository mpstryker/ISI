global Degperpix xwinpix ywinpix azcp elcp

function [xdeg ydeg] = p2d(xpix, ypix)
    xdeg = ((xpix-(xwinpix/2))/degperpix)+cpaz;
    ydeg = ((ypix-(ywinpix/2))/degperpix)+cpel;
end
