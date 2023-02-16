function Pfeil_1
global n_elems k my sigma masses dt th M elem_pos_x speed rect_len reset_loop

p = 15;
g = 0.05;
q = 1;
th = 0.2;

om=1;
a=0.063162931;

T = 200;
M = 500;

dt = T / (M-1);

n_elems = 25;
len = 10;

k = 0.5;
my = -0.05;
sigma = 100;
speed = 0.5;


reset_loop = 0;



% empiristic model for flight path in w direction
th_f = @(t) exp(-t(:)/10).*sin(t(:));
n0 = @(x) sin(pi*x(:)')/10-a;
w = @(x,t) sin(om*t(:))*th_f(t(:))*n0(x(:)')+p+q*t(:)-g*t(:).^2/2+0.05*(q/g-t(:)')*x(:)';


x=linspace(0,len,n_elems);
h=x(2)-x(1);
xc=x(1:n_elems)+h/2;



% start conditions
masses = ones(n_elems, 1);

elem_pos_x = linspace(0, len, n_elems)';
elem_vel_x = zeros(n_elems, 1);

elem_pos_y = linspace(0, len, n_elems)';
elem_pos_y = ones(n_elems, 1) * 0;
elem_vel_y = zeros(n_elems, 1);
elem_vel_y(round(n_elems/2)) = -6;


elem_vel_y = elem_vel_y - mean(elem_vel_y);
%elem_vel_y = sin(2*pi*elem_pos_x/(len*2))*0.5;

X0 = cat(1, elem_pos_x, elem_vel_x, elem_pos_y, elem_vel_y);
tspan = [0,T];


% Theta method
t = linspace(0, T, M);
Y = Theta(X0)';
% ode method
% [t,Y] = ode45(@RechteSeite, tspan, X0);


rect_len = (len/n_elems)/2;

% figure preparation
h1 = figure(1); close(h1); h1 = figure(1);
h1.WindowState = 'maximized';

ax = subplot(1,3,3);
ax.XLim = [0 100];
ax.YLim = [-50 50];
ax.ZLim = [0 40];

ax = subplot(1,3,1);
hold on
for i = 1:n_elems
    plot(ax, t, Y(:,n_elems + i))
end
hold off

set(h1,'KeyPressFcn',@key_detect_func);


while true
    for r = 1:size(t,2)
        if(not(ishandle(h1)))
            close all
            return;
        else
        if(reset_loop == 1)
            reset_loop = 0;
            break;
        end
        sgtitle(sprintf("%d/%d",r,size(t,2)));
        view(2)
        ax = subplot(1,3,2);
        cla
        hold on
        average = mean(Y(r,1:n_elems));
        for i = 1:n_elems
            x = Y(:,i);
            y = Y(:,n_elems+i);
            rectangle(ax, "Position", [x(r)-rect_len/2, y(r)-rect_len/2, rect_len, rect_len])
        end
        text(10,10,sprintf("%d/%d",r,size(t,2)))
        ax = h1.CurrentAxes;
        ax.XLim = [average-7 average+7];
        ax.YLim = [-5 7];
        
        radius = 0.1;
        c = [10.5 -0.5];
        pos = [c-radius 2*radius 2*radius];
        rectangle('Position',pos,'Curvature',[1 1])
        
        %
        %size(Y(r,1:n_elems)')
        %size(Y(r,n_elems+1:2*n_elems)')
        %size(w(t(k)*ones(1,n_elems)))
        w0=w(xc,t(r)/5);
        
        hold off
        ax = subplot(1,3,3);
        cla
        hold on
        view([r/20-20 90-r/5])
        grid on
        plot3(ax, Y(r,1:n_elems)', Y(r,n_elems+1:2*n_elems)', w0);
        [cy_X,cy_Y,cy_Z] = cylinder(0.1);
        cy_Z = cy_Z*20;
        cy_X = cy_X + 10.5;
        cy_Y = cy_Y - 0.5;
        surf(ax, cy_X, cy_Y, cy_Z);
        
        %ax.XLim = [0 100];
        %ax.YLim = [-50 50];
        %ax.ZLim = [0 40];
        drawnow;
        hold off
        
        pause(0.0001)
        end
    end
end
end

function key_detect_func(src, event)
global reset_loop
if(event.Key == 114)
    reset_loop = 1;
end
end


function vsav = Theta(VV0)
global n_elems dt th k my sigma masses M elem_pos_x speed rect_len

%kappa = masses * 1/k;
kappa = k./masses;

Z = zeros(n_elems);
I = diag(ones(n_elems,1));

L = diag(ones(n_elems,1)*-2) + diag(ones(n_elems-1, 1),1) + diag(ones(n_elems-1, 1),-1);
L(1,1) = -1;
L(n_elems,n_elems) = -1;

L = kappa .* L;

D = (my./masses) .* I;

e = ones(n_elems, 1);
B = full(spdiags([-1*e 4*e -6*e 4*e -1*e], -2:2, n_elems, n_elems));
B(1:2,1:2) = [-1,2;2,-5];
B(end-1:end,end-1:end) = [-5,2;2,-1];

B = B * sigma;

A = [Z,Z,Z,Z;
    Z,Z,Z,Z;
    Z,Z,Z,I;
    Z,Z,(L+B)./masses,D];


VV = VV0;
v = VV(1:n_elems*2);
vsav = v;
AL = speye(2*2*n_elems) - (1-th)*dt*A;
AR = speye(2*2*n_elems) + th*dt*A;

z = zeros(n_elems,1);

F = [ones(n_elems,1)*speed;z;z;z];

for k=2:M
    VV = AL \ (AR*VV + dt*F);
    v = cat(1, VV(1:n_elems), VV(n_elems*2+1:n_elems*3));
    ax = gca;
    cla(ax);
    vsav=[vsav,v];
end
end




function dY = RechteSeite(t,X)
global n_elems k my sigma masses speed rect_len

kappa = masses * 1/k;


Z = zeros(n_elems);
I = diag(ones(n_elems,1));

L = diag(ones(n_elems,1)*-2) + diag(ones(n_elems-1, 1),1) + diag(ones(n_elems-1, 1),-1);
L(1,1) = -1;
L(n_elems,n_elems) = -1;

L = kappa .* L;

D = (my./masses) .* I;

e = ones(n_elems, 1);
B = full(spdiags([-1*e 4*e -6*e 4*e -1*e], -2:2, n_elems, n_elems));
B(1:2,1:2) = [-1,2;2,-5];
B(end-1:end,end-1:end) = [-5,2;2,-1];

B = B * sigma;
B = B./rect_len^4;

A = [Z,Z,Z,Z;Z,Z,Z,I;Z,Z,Z,Z;Z,(L+B)./masses,Z,D];

%A = [ZZ, ZZ+I; ZZ+(L+B)./masses, ZZ+D];
b = zeros(n_elems*4,1);
%b(1:n_elems,1) = 1;
dY = A*X + b;
end


