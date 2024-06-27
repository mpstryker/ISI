function udpmon
udp = pnet('udpsocket',8936);

%data=pnet(con,'read' [,size] [,datatype] [,swapping] [,'view'] [,'noblock'])
done = 0;

while ~done
    
    size=pnet(udp,'readpacket'); % size returns 8 (bytes)

    %data=pnet(udp,'read' , 6, 'byte');
    data = zeros(1:8);
    data=pnet(udp,'read', size);
    fprintf(1, '%c %c %x %x %x %x %x %x \n', data(1), data(2), data(3), data(4), data(5), data(6), data(7), data(8) );

    [keyIsDown, secs, keyCode] = KbCheck;
    if keyIsDown
        done=1;
    end
end


pnet('closeall')