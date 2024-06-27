function xdeg = p2dx(xpix)
global degperpix xwinpix azcp;
    xdeg = ((xpix-(xwinpix/2))/degperpix)+ azcp;
end
