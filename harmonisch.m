function harmonisch
    X0 = [-3;-1]; tspan = [0,15];
    [t,X] = ode45(@RechteSeite,tspan,X0);
    w = X(:,1);
    wexakt = -sin(t)-cos(t)-2;
    plot(t,w,'r*',t,wexakt,'b')
end

function dX = RechteSeite(t,X)
    A = [0,1;-1,0]; b = [0;-2];
    dX = A*X+b;
end