function harmonisch_n_massenschwinger_quer
    global n_elems
    n_elems = 50;
    len = 10;
    
    elem_pos_y = linspace(0, len, n_elems)';
    elem_pos_x = linspace(0, len, n_elems);

    elem_pos_y = ones(n_elems, 1) * 0;
    
    elem_vel_y = zeros(n_elems, 1);

    elem_vel_y = sin(2*pi*elem_pos_x'/(len*2))*0.25;
    
    X0 = cat(1, elem_pos_y, elem_vel_y);
    
    middle_vel_elem = n_elems+round(n_elems/2);
    
    
    %X0(middle_vel_elem-1) = -0.5;
    %X0(middle_vel_elem) = -1;
    %X0(middle_vel_elem+1) = -0.5;

    
    tspan = [0,500];
    [t,Y] = ode45(@RechteSeite, tspan, X0);
   
    figure(1);
    cla
    hold on
    for i = 1:n_elems
        plot(t, Y(:,i))
    end
    hold off
       
    h1 = figure(2);
    
    rect_len = (len/n_elems)/2;

    for r = 1:size(t)
        if(not(ishandle(h1)))
            return
        end
        ax = h1.CurrentAxes;
        cla(ax);
        for i = 1:n_elems
            y = Y(:,i);
            rectangle(ax, "Position", [elem_pos_x(i)-rect_len/2, y(r)-rect_len/2, rect_len, rect_len])
        end
        ax = h1.CurrentAxes;
        ax.XLim = [-1 11];
        ax.YLim = [-5 7];

        pause(0.5/n_elems)
    end
    
end



function dY = RechteSeite(t,X)
    global n_elems

    k = 1;
    my = -0.05;
    sigma = 16;
    
    masses = ones(n_elems, 1);
    %masses(2) = 20;
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
    
    A = [Z, I; (L+B)./masses, D];
   
    dY = A*X;
end

