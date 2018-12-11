% Solves 2D Heat eqn
%     T_t = alpha*(T_xx +T_yy)
% with using ADI scheme.
%%
clearvars;close all;clc;
tic
%% parameters
lx=3.5;                     % length 
ly=3.5;             
alpha=0.654;                % diffusion constant

dx=0.1;dy=dx;               % space step
mx=lx/dx;                   % number
my=ly/dy;
nx=mx+1;
ny=my+1;
x=linspace(0,lx,nx);
y=linspace(0,ly,ny);

dt=0.001;                   % time step
tmax=0.5;                   % final time
beta2=(dx*dx)/(dy*dy);      % constants 
r=alpha*dt/(dx*dx);

t=0;                        % initialize
k=1;
T=zeros(nx,ny);             
Thalf=zeros(nx,ny);
%% BCs

%left 
T(1,:)=200;
% right
T(nx,:)=0;
% bottom
T(:,1)=200;
% top 
T(:,ny)=0;

% set up TDMs
a1(2:mx)=1;
b1(1:mx+1)=-2-2/r;
c1(1:mx+1)=1;
d1(1:mx+1)=0;
a1(1)=0;b1(1)=1;c1(1)=0;
a1(mx+1)=0;b1(mx+1)=1;c1(mx+1)=0;
a2(1:my+1)=beta2;

b2(1:my+1)=-2*beta2-2/r;
c2(1:my+1)=beta2;
a2(1)=0;b2(1)=1;c2(1)=0;
a2(my+1)=0;b2(my+1)=1;c2(my+1)=0;
d2(1:my+1)=0;


% UU={};

%%
while (t<tmax)
    for j=2:my
        % call TDMA
        d1(1)=T(1,j);
        % Dirichlet
        d1(mx+1)=T(mx+1,j);
        %Neumann
        d1(mx+1)=T(mx+1,j)+ r*(T(mx,j)-2*T(mx+1,j)+T(mx,j));
        for i=2:mx
            d1(i)=-beta2*T(i,j-1)-beta2*T(i,j+1)+(2*beta2-2/r)*T(i,j);
        end
        Thalf(:,j)=solver_tdma(mx+1,a1,b1,c1,d1);
    end
    for i=2:mx
        % call TDMA
        d2(1)=T(i,1);
        d2(my+1)=T(i,my+1);
        for j=2:my
            d2(j)=-Thalf(i-1,j)-Thalf(i+1,j)+(2-2/r)*Thalf(i,j);
        end
        T(i,:)=solver_tdma(my+1,a2,b2,c2,d2);
    end
    t=t+dt;
%     UU{k}=T';

    if mod(k,50)==0
        h=contour(x,y,T,'Fill','on');colorbar;
        title(sprintf('time=%0.5f hr',t));
        drawnow
    end
    k=k+1;
    
    
end
elapsedTime=toc;
fprintf('Took %0.4f seconds \n',elapsedTime);
%% plots
figure(1)
subplot(2,1,1)
contour(x,y,T)
title(sprintf('Temperature at t=%0.2f hr',tmax));
xlabel('x');ylabel('y');colorbar;
subplot(2,1,2)
pcolor(x,y,T);shading interp,
title(sprintf('Temperature at t=%0.2f hr',tmax));
xlabel('x');ylabel('y'),colorbar;

%%
function x = solver_tdma(n,a,b,c,f)
% Solve an nxn tridiagonal system with sub-diagonal a, diagonal b,
% superdiagonal c, and rhs f
x = zeros(size(f));
c(1)=c(1)/b(1);
f(1)=f(1)/b(1);
% Forward elimination
for i=2:n
    p = 1.0/(b(i)-c(i-1)*a(i));
    c(i) = c(i)*p;
    f(i) = (f(i) - a(i)*f(i-1))*p ;
end
% Back substitution
x(n) = f(n);
for i=n-1:-1:1
    x(i) = (f(i) - c(i)*x(i+1));
end
end
