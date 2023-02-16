function KetteLaeng
    global AA bb
    k = 1; m = 1; r = 0; L = 9;
    A = (k/m)*[-2,1;1,-2]; Z = [0,0;0,0]; I = [1,0;0,1]; AA = [Z,I;A,-r*I];
    b = (k/m)*[0;L]; z = [0;0]; bb = [z;b];
    u1_0 = 0.75; u2_0 = 2.25; u1s_0 = 0; u2s_0 = 0;
    U0 = [u1_0;u2_0]; U1 = [u1s_0;u2s_0]; UU0 = [U0;U1];
    tspan = [0,10];
    [t,UU] = ode45(@RechteSeite,tspan,UU0);
    u1     = UU(:,1);
    u2     = UU(:,2);
    h1 = figure(1); close(h1); h1 = figure(1);
    set(h1,'Position',[10 100 500 500]);
    plot(t,u1,'r',t,u2,'b',t,ones(size(t)),'k:',t,2*ones(size(t)),'k:','LineWidth',3)
    ylim([0,L]);
    legend('u1','u2','Location','best')
    xlabel('t')
    ylabel('u1,u2')
    title('Kette: laengliche Auslenkungen')
    text(2,2.8,'ohne Dämpfung','FontName','Times','FontSize',15,'FontWeight','bold')
    set(gca,'FontName','Times','FontSize',15,'FontWeight','bold')
end

function F = RechteSeite(t,UU)
    global AA bb
    F = AA*UU + bb;
end
