% Solves advection equation
%     u_t+c*u_x = 0
% with 1st O Upwind and Lax-Wendroff Schemes
%%
% Parameters
courants=[0.45 0.9 1 2];            % courant numbers
dx=5;                               % space step
c=400;                              % propagation speed
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

for kk=1:4
    
    co=courants(kk);                % courant number co=(c*dt)/dx
    dt=(co*dx)/c;                   % time step
    t=0:dt:0.45;                    % time range
    nt=length(t);                   % number of time steps
    
    %initialize
    uupwind=zeros(nt,nx); uexact=zeros(nt,nx);
    uupwind(1,:)= u0;  uexact(1,:)=u0;
    ulax=uupwind;
    
    for j=2:nt
        %bc
        uupwind(j,1)=0;uupwind(j,end)=0;
        uexact(j,1)=0;uexact(j,end)=0; 
        
        for i=2:nx-1
            %Upwind
            uupwind(j,i)=uupwind(j-1,i)-co*(uupwind(j-1,i)-uupwind(j-1,i-1));
            %Lax Wendroff
            ulax(j,i)=ulax(j-1,i)-(co/2)*(ulax(j-1,i+1)-ulax(j-1,i-1))+(co^2/2)*(ulax(j-1,i+1)-2*ulax(j-1,i)+ulax(j-1,i-1));
        end 
        
        %exact
        for v=1:nx
            if((x(v)-c*t(j))>=50 && (x(v)-c*t(j))<110)
                uexact(j,v)=100*sin(pi*((x(v)-c*t(j)-50)/60));
            end
        end
    end
    
    figure(1);
    hold on
    plot(x,uupwind(end,:));
    
    figure(2);
    hold on
    plot(x,ulax(end,:));
    
end
%%
% plot
hold on
figure(1);
plot(x,uexact(end,:),'--'); 
legend('Courant Number= 0.45','Courant Number= 0.9','Courant Number= 1','Courant Number= 2','Exact','Location','Northwest');
title('1st order Upwind scheme')
xlabel('x');ylabel('u');
axis([0 400 0 100]);

figure(2);
plot(x,uexact(end,:),'--'); 
legend('Courant Number= 0.45','Courant Number= 0.9','Courant Number= 1','Courant Number= 2','Exact','Location','Northwest');
title('Lax-Wendroff scheme')
xlabel('x');ylabel('u');
axis([0 400 -15 105]);
