global Degperpix xwinpix ywinpix azcp elcp

function xpix = d2px(xdeg)
    xpix = ((xdeg-cpaz)*degperpix)-(xwinxpix/2);
end
