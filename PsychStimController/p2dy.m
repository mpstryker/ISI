function ydeg = p2dy(ypix)
global ywinpix degperpix elcp;
    ydeg = ((ypix-(ywinpix/2))/degperpix)+elcp;
end

