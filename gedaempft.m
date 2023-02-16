function Gedaempft
X0 = [-3;-1]; tspan = [0,15];
[t,X] = ode45(@RechteSeite,tspan,X0);
w = X(:,1);
v = X(:,2);
wexakt = -exp(-t/2).* (cos(sqrt(3)*t/2)+sqrt(3)*sin(sqrt(3)*t/2))-2;
plot(t,w,'r*',t,wexakt,'b', t, v, 'g*')
end
function dX = RechteSeite(t,X)
A = [0,1;-1,-1]; b = [0;-2];
dX = A*X+b;
end