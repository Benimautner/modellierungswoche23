function harmonisch_zweimassenschwinger_2
    X0 = [3.7;6.1;0;0]; 
    tspan = [0,50];
    [t,X] = ode45(@RechteSeite, tspan, X0);
    
    x1 = X(:,1);
    x2 = X(:,2);
    figure(1);
    plot(t, x1, "r", t, x2, "g");
    
    h1 = figure(2);
    
    rect_len = 1;
    
    limx = [min(cat(1,x1,x2))-rect_len * 2 max(cat(1,x1,x2))+rect_len * 2];
    
    
    for r = 1:size(t)
        if(not(ishandle(h1)))
            return
        end
        ax = h1.CurrentAxes;
        cla(ax);
        rectangle(ax, "Position", [x1(r)-rect_len/2, -rect_len/2, rect_len, rect_len])
        rectangle(ax, "Position", [x2(r)-rect_len/2, -rect_len/2, rect_len, rect_len])
        ax = h1.CurrentAxes;
        ax.XLim = limx;
        range = ax.XLim(2) - ax.XLim(1);
        ax.YLim = [-range/2 range/2];
        pause(0.1)
    end
    
end

function dX = RechteSeite(t,X)
    k1 = 1;
    k2 = 1;
    k3 = 1;
    m1 = 1;
    m2 = 1;
    u1 = 3;
    u2 = 6;
    my = 0;
    
    A = [0,           0,            1,   0;
         0,           0,            0,   1;
         (-k1-k2)/m1, k2/m1,       -my,  0;
         k2/m2,       (-k2-k3)/m2,  0,  -my];
    b = [0; 0; (k1*u1+k2*u1-k2*u2)/m1; (k2*u2-k2*u1+k3*u2)/m2];
    dX = A*X + b;
end