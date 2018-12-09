%%
% Parameters
co=0.7;                    % courant number co=(c*dt)/dx
dx=5;                      % space step
x=0:dx:400;                % x range
nx=length(x);              % number of space steps
%%
% ic
u0=zeros(1,nx);
for k=1:nx
    %propagation speed c
    if(x(k)>=0 && x(k)<=200 )
        c=400;
    else
        c=-400;
    end
    %initial condition
    if((x(k)>=50 && x(k)<110) || (x(k)>=290 && x(k)<350) )
        u0(k)=100*sin(pi*((x(k)-50)/60));
    end
end

dt=abs((co*dx)/c);         % time step
t=0:dt:0.3;                % time range
nt=length(t);              % number of time steps

%initialize
ulax=zeros(nt,nx); uexact=zeros(nt,nx);
ulax(1,:)=u0; uexact(1,:)=u0;
uupwind=ulax;
%%

for j=2:nt
    %boundary conditions
    ulax(j,1)=0;ulax(j,end)=0;
    uupwind(j,1)=0;uupwind(j,end)=0;
    uexact(j,1)=0;uexact(j,end)=0;
    
    for i=2:nx-1
        %Lax-Wendroff
        if(x(i)>=0 && x(i)<=200 )
            c=400;
        else
            c=-400;
        end
        co=c*(dt/dx);
        
        ulax(j,i)=ulax(j-1,i)-(co/2)*(ulax(j-1,i+1)-ulax(j-1,i-1))+(co^2/2)*(ulax(j-1,i+1)-2*ulax(j-1,i)+ulax(j-1,i-1));
        
        %Upwind
        if(x(i)>=0 && x(i)<=200 )
            uupwind(j,i)=uupwind(j-1,i)-co*(uupwind(j-1,i)-uupwind(j-1,i-1));
        else
            uupwind(j,i)=uupwind(j-1,i)-co*(uupwind(j-1,i+1)-uupwind(j-1,i));
        end 
    end
    %exact
    for v=1:nx
        if(((x(v)-c*t(j))>=50 && (x(v)-c*t(j))<110) || ((x(v)-c*t(j))>=290 && (x(v)-c*t(j))<350) )
            uexact(j,v)=100*sin(pi*((x(v)-c*t(j)-50)/60));
        end
    end  
end
%%
figure(1)
plot(x,ulax(j,:),'LineWidth',.5);
hold on
plot(x,uupwind(end,:),'LineWidth',.5);
plot(x,uexact(end,:),'g','LineWidth',.5);
% axis([0 400 0 100]);
title('u(x) at t=0.30s for Lax-Wendroff, Upwind, and exact solution');
xlabel('x');ylabel('u');
legend('Lax-Wendroff','Upwind','Exact');