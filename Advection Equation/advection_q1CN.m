% Solves advection equation
%     u_t+c*u_x = 0
% with Crank-Nicolson scheme using either 1st O BD or 2nd O CD for space
%%
clearvars;clc;
hold on;
selection=input(' Press 1 - Crank-Nicolson 1st O BD,\n Press 2 - Crank-Nicolson 2nd O CD\n');
%%
%Parameters
courants=[0.45 0.9 1 2];            % courant numbers
dx=5;                               % space step
c1=400;                             % propagation speed
x=0:dx:400;                         % x range
nx=length(x);                       % number of space steps
%%
% Initial condition
u0=zeros(1,nx);
for k=1:nx
    if(x(k)>=50 && x(k)<110)
        u0(k)=100*sin(pi*((x(k)-50)/60));
    end
end
%%
% loop courant numbers
for kk=1:4
    co=courants(kk);                % courant number co=(c*dt)/dx
    dt=(co*dx)/c1;                  % time step
    t=0;tmax=0.45;                  % time
    
    %initial condition
    u=u0;
    
    % tridiagonal components
    a(1:nx)=-0.25*co;
    b(1:nx)=1;
    c(1:nx)=0.25*co;
    
    % boundary conditions
    a(1)=0;b(1)=1;c(1)=0; 
    a(nx)=0;b(nx)=1;c(nx)=0;
    d(1:nx)=0;
    
    %time loop
    while t<tmax
        for j=2:nx-1
            if(selection==2)
                %2nd Order CD
                d(j)=u(j)-0.25*co*(u(j+1)-u(j-1));
            else
                %1st Order BD
                d(j)=u(j)-0.5*co*(u(j)-u(j-1));
            end
        end
        d(1)=0;d(nx)=0;
        u = solver_tdma(nx,a,b,c,d);

        t=t+dt;
    end   
    plot(x,u);
end
%%
% Exact Solution
uexact=zeros(1,nx);
for v=1:nx
    if((x(v)-c1*0.45)>=50 && (x(v)-c1*0.45)<110)
        uexact(v)=100*sin(pi*((x(v)-c1*0.45-50)/60));
    end
end
%%
% plot
hold on
plot(x,uexact);
if(selection==1)
    axis([0 400 0 100]);
    title('Crank Nicolson 1st Order Backward Differnce')
else
    axis([0 400 -20 110]);
    title('Crank Nicolson 2nd Order Central Difference')
end
legend('Courant Number= 0.45','Courant Number= 0.9','Courant Number= 1','Courant Number= 2','Exact','Location','Northwest');
xlabel('x');ylabel('u');
%%
%TDMA solver
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
