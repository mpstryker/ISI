global Degperpix xwinpix ywinpix azcp elcp

function [xpix ypix] = d2p(xdeg, ydeg)
    xpix = ((xdeg-cpaz)*degperpix)-(xwinxpix/2);
    ypix = ((ydeg-cpel)*degperpix)-(ywinypix/2);
end

