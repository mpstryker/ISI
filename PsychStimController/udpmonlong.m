function udpmon
udp = pnet('udpsocket',8936);

%data=pnet(con,'read' [,size] [,datatype] [,swapping] [,'view'] [,'noblock'])
done = 0;

while ~done
    
    size=pnet(udp,'readpacket'); % size returns 8 (bytes)

    %data=pnet(udp,'read' , 6, 'byte');
    data = zeros(1:2,'uint32');
    data=pnet(udp,'read', size/4, 'uint32');
    fprintf(1, '%u %u \n', data(1), data(2) );
    

    [keyIsDown, secs, keyCode] = KbCheck;
    if keyIsDown
        done=1;
    end
end


pnet('closeall')