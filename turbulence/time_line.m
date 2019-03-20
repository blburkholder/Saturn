function [] = time_line(rbd,rx,i)
    if rbd(7,i) == 2
        h = patch([rx(i),rx(i),rx(i+1),rx(i+1)],[0.5,0,0,0.5],'g');
        set(h,'EdgeColor','none')
    elseif rbd(7,i) == -1
        h = patch([rx(i),rx(i),rx(i+1),rx(i+1)],[0,-0.5,-0.5,0],'b');
        set(h,'EdgeColor','none')
    elseif rbd(7,i) == 1
        h = patch([rx(i),rx(i),rx(i+1),rx(i+1)],[0.5,1.0,1.0,0.5],'r');
        set(h,'EdgeColor','none')
    elseif rbd(7,i) == -2
        h = patch([rx(i),rx(i),rx(i+1),rx(i+1)],[-0.5,-0.75,-0.75,-0.5],'k');
        set(h,'EdgeColor','none')
    end