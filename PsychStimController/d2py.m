global Degperpix xwinpix ywinpix azcp elcp

function ypix = d2py(ydeg)
    ypix = ((ydeg-cpel)*degperpix)-(ywinypix/2);
end
