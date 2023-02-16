function pfeil
    global n_elems k my sigma masses dt th M elem_pos_x
    
    p = 100;
    g = 0.05;
    q = 1;
    th = 0.2;
    
    T = 400;
    M = 500;
    
    dt = T / (M-1);
    
    u = @(x,t) x(:) - q*t(:);
    w = @(t) p - g*t(:).^2/2;
       
    n_elems = 50;
    len = 10;
    
    k = 1;
    my = -0.05;
    sigma = 2;
    
    masses = ones(n_elems, 1);
    
    elem_pos_y = linspace(0, len, n_elems)';
    elem_pos_x = linspace(0, len, n_elems)';
    
    

    elem_pos_y = ones(n_elems, 1) * 0;
    
    elem_vel_x = zeros(n_elems, 1);
    elem_vel_y = zeros(n_elems, 1);
    
    elem_vel_y(round(n_elems/2)) = -8;
    a = mean(elem_vel_y);
    elem_vel_y = elem_vel_y - a;

    %elem_vel_y = sin(2*pi*elem_pos_x/(len*2))*0.5;
    
    X0 = cat(1, elem_pos_x, elem_vel_x, elem_pos_y, elem_vel_y);
    
    %middle_vel_elem = n_elems+round(n_elems/2);
    
    tspan = [0,T];
    t = linspace(0, T, M);
    Y = Pheta(X0)';
    %[t,Y] = ode45(@RechteSeite, tspan, X0);
    
    size(Y)
    size(t)
   
    figure(1);
    cla
    hold on
    for i = 1:n_elems
        plot(t, Y(:,n_elems + i))
    end
    hold off
       
    h1 = figure(2);
    
    rect_len = (len/n_elems)/2;
    
    results = Y(:,1:n_elems);
    
    mode = 0;
    
    hold on
    size(t)
    while true
        for r = 1:size(t,2)
            if(not(ishandle(h1)))
                return
            end
            if(mode == 0)
                ax = h1.CurrentAxes;
                cla(ax);
                for i = 1:n_elems
                    x = Y(:,i);
                    y = Y(:,n_elems+i);
                    rectangle(ax, "Position", [x(r)-rect_len/2, y(r)-rect_len/2, rect_len, rect_len])
                end
                text(10,10,sprintf("%d/%d",r,size(t,2)))
                ax = h1.CurrentAxes;
                ax.XLim = [-1 11];
                ax.YLim = [-5 7];
            else
                %
                cla
                %size(Y(r,1:n_elems)')
                %size(Y(r,n_elems+1:2*n_elems)')
                %size(w(t(k)*ones(1,n_elems)))


                plot3(Y(r,1:n_elems)', Y(r,n_elems+1:2*n_elems)', w(t(k)*ones(1,n_elems)));
    %             ax = gca;
    %             ax.XLim = [-10 20];
    %             ax.YLim = [-5 5];
    %             ax.ZLim = [99 101];
                drawnow;
            end
            pause(0.3/size(t,2))
        end
    end
    hold off
    
    
end


function vsav = Pheta(VV0)
    global n_elems dt th k my sigma masses M elem_pos_x

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
    
    F = [ones(n_elems,1);z;z;z];
    
    for k=2:M
        VV = AL \ (AR*VV + dt*F);
        v = cat(1, VV(1:n_elems), VV(n_elems*2+1:n_elems*3));
        plot(elem_pos_x, v(n_elems-1:n_elems-1)); drawnow;
        vsav=[vsav,v];
    end
end




function dY = RechteSeite(t,X)
    global n_elems k my sigma masses

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
    
    A = [Z,Z,Z,Z;Z,Z,Z,I;Z,Z,Z,Z;Z,(L+B)./masses,Z,D];
    
    %A = [ZZ, ZZ+I; ZZ+(L+B)./masses, ZZ+D];
    b = zeros(n_elems*4,1);
    %b(1:n_elems,1) = 1;
    dY = A*X + b;
end


