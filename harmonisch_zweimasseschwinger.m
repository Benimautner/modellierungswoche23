function harmonisch_zweimasseschwinger
    X0 = [-1;1;0;0]; 
    tspan = [0,15];
    [t,X] = ode45(@RechteSeite, tspan, X0);
    
    x1 = X(:,1);
    x2 = X(:,2);
    
    plot(t, x1, "r*", t, x2, "g*");
    
    h1 = figure(1); close(h1); h1 = figure(1);
    
    axis([-2 2 -2 2])
    rect_len = 0.1
    for r = 1:size(t)
        cla
        rectangle("Position", [x1(r)-rect_len/2, -rect_len/2, rect_len, rect_len])
        rectangle("Position", [x2(r)-rect_len/2, -rect_len/2, rect_len, rect_len])
        pause(0.1)
    end
    
end

function dX = RechteSeite(t,X)
    c1 = 1;
    c2 = 1;
    c3 = 1;
    m1 = 1;
    m2 = 1;
    A = [0,0,1,0;
         0,0,0,1;
         ((-c1-c2)/m1),(c2/m1),0,0;
         (c2/m2),((-c2-c3)/m2),0,0];
    dX = A*X;
end