N=50;
M=50;
T=12.79;

a=0.063162931;p=0;q=64;g=10;om=1;
if(p+q*T-g*T^2/2<0)
    error('aufpralhöhe negativ!');
end
n0=@(x)sin(pi*x(:)')/10-a;
n1=@(x)pi*cos(pi*x(:)')/10;
th=@(t)exp(-t(:)/10).*sin(t(:));
u=@(x,t)x(:)'+ q*t(:);
w=@(x,t)sin(om*t(:))*th(t(:))*n0(x(:)')+p+q*t(:)-g*t(:).^2/2+0.05*(q/g-t(:)')*x(:)';
v=@(x,t)cos(om*t(:))*th(t(:))*n0(x(:)');
vx=@(x,t)cos(om*t(:))*th(t(:))*n0(x(:)');

x=linspace(0,1,N);h=x(2)-x(1);
xc=x(1:N-1)+h/2;
t=linspace(0,T,M)';

yc=[];
for k=1:M
   
    u0=u(xc,t(k));
    v0=v(xc,t(k));
    v1=vx(xc,t(k));
    w0=w(xc,t(k));
    a0=th(t(k))*n0(xc);
    a1=th(t(k))*n1(xc);
    L=h*sum(sqrt(1+a1.^2));
    X=h*sum(xc.*sqrt(1+a1.^2))/L;
    Y=h*sum(a0.*sqrt(1+a1.^2))/L;
    yc=[yc,Y];
    
     subplot(1,3,1)
    %plot(xc,a0,'b',xc,v0,'m',xc,w0,'c:',X,Y,'r*');ylim([-0.5,+0.5]);xlim([0,1]);
    xlabel('x');ylabel('v')
    title(['t=',num2str(t(k))])
    set(gca,'FontName','Times','FontSize',15,'FontWeight','bold');
    pbaspect([1 1 1])
    drawnow;
    
    subplot(1,3,2)
    plot3(u0,v0,w0,'b',...
        u(0.5,t(k)),Y,w(0.5,t(k)),'r*');
    xlabel('u');ylabel('v');zlabel('w')
    umin=min(u(xc,t(k)));umax=max(u(xc,t(k)));
    wmin=min(w(xc,t(k)));wmax=max(w(xc,t(k)));
    axis([umin-1,umax+1,-1.5,1.5,wmin-1,wmax+1])
    grid on;
    view ([k 30])
    set(gca,'FontName','Times','FontSize',15,'Fontweight','bold');
    title(['t=',num2str(t(k))])
    pbaspect([1 1 1])
    drawnow;
end
subplot(1,3,3)
size(t(size(yc)))
plot(t(1:length(yc)),yc,'r',t,0*zeros(size(t)),'k:')
xlabel('t')
ylabel('vbar')
title('Schwerpunkt')
set(gca,'Fontname','Times','FontSize',15,'FontWeight','bold');
pbaspect([1 1 1])