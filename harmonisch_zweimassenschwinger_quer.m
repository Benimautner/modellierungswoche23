function harmonisch_zweimassenschwinger_quer
    X0 = [0.1;-0.1;0;0]; 
    tspan = [0,20];
    [t,Y] = ode45(@RechteSeite, tspan, X0);
    
    y1 = Y(:,1);
    y2 = Y(:,2);
    figure(1);
    plot(t, y1, "r", t, y2, "g");
    
    h1 = figure(2);
    

    rect_len = 1;

    limx = [min(cat(1,y1,y2))-rect_len * 2 max(cat(1,y1,y2))+rect_len * 2];
    
    for r = 1:size(t)
        if(not(ishandle(h1)))
            return
        end
        ax = h1.CurrentAxes;
        cla(ax);
        rectangle(ax, "Position", [3, y1(r)-rect_len/2, rect_len, rect_len])
        rectangle(ax, "Position", [6, y2(r)-rect_len/2, rect_len, rect_len])
        ax = h1.CurrentAxes;
        ax.XLim = [0 10];
        range = ax.XLim(2) - ax.XLim(1);
        ax.YLim = limx;
        pause(0.1)
    end
    
end

function dX = RechteSeite(t,X)
    k1 = 1;
    k2 = 1;
    k3 = 1;
    m1 = 1;
    m2 = 1;
    l = 3;
    my = 0.5;
    A = [0,0,1,0;
         0,0,0,1;
         (-k1*3)/(2*l*m1)-(k2*3)/(2*l*m1), (k2*3)/(2*l*m1), -my/m1, 0;
         (k2*3)/(2*l*m2), -(k3*3)/(2*l*m2)-(k2*3)/(2*l*m2), 0, -my/m2];

    dX = A*X;
end

