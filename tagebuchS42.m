p = 10; q = 1; g = 0.1; u = @(x,t) x(:)' + q*t(:); w = @(t) p - g*t(:).^2/2;
N = 50; M = 100; T = 20; dt = T/(M-1); dx = 1/(N-1);
x = linspace(0,1,N)'; t = linspace(0,T,M);
sg = 1.0e-2; ka = 1.0e-2; mu = 1.0e-1; th = 0.5;
m = ones(N,1); M1 = spdiags(1./m,0,N,N); e = ones(N,1);
B = spdiags([-e 4*e -6*e 4*e -e],-2:2,N,N);
B(1,1) =-1; B(1,2) = 2; B(2,1) = 2; B(2,2) =-5;
B(N-1,N)= 2; B(N-1,N-1)=-5; B(N,N) =-1; B(N,N-1)= 2; B=sg*B/dx^4;
L = spdiags([e -2*e e],-1:1,N,N); L(1,1)=-1; L(N,N)=-1; L=ka*L/dx^2;
Z = 0*speye(N); I = speye(N); AA = [Z,I;M1*(B+L),-mu*M1];
V0 = zeros(N,1); V1 = zeros(N,1);
z = zeros(N,1); F = [z;z];
b = 0.5; a = 1.96*b/pi;
V1 = a - b*sin(pi*x(:));
VV0 = [V0;V1];
VV = VV0;
AL = speye(2*N) - (1-th)*dt*AA;
AR = speye(2*N) + th*dt*AA;
v = VV(1:N); vsav = v;
subplot(1,2,1)
plot(x,v); drawnow;
for k=2:M
VV = AL \ (AR*VV + dt*F);
v = VV(1:N); vsav = [vsav,v];
subplot(1,2,1)
hold on; plot(x,v);
drawnow; hold off;
subplot(1,2,2)
plot3(u(x,t(k)),v,w(t(k))*ones(1,N));
grid on; view([k 30]); drawnow;
end