% Solves Lid Driven cavity flow 
% using streamfunction-vorticity formulation
% with forward time explicit scheme
% given streamfunction boundary conditions 
% left Inlet and right outlet.
% Implements upwind blending.
%%
clearvars;close all;clc;tic
%% parameters
h=.01;              % space step (dx=dy=h)
L=0.3;              % Length
n=L/h+1;            % number
x=0:h:L;
y=0:h:L;
U=5;                % Lid Velocity
dt=0.002;           % time step
endTime=3;          % final time

Re=200;             % Reynolds no.
nu=(U*L)/Re;        % Courant no.

vort=zeros(n,n);    % initialize
psi=zeros(n,n);
s=0;
t=0;
iter=0;
q=0.5;              % blending factor
%% initial
for i=1:n
    for j=1:n
        psi(i,j,1)=0;
        vort(i,j,1)=0;
        if j==n
            psi(i,j,1)=0;
            vort(i,j,1)=-2*U/h;
        end
    end
end

%% time loop
while t<endTime
    
    vortold=vort;
    psiold=psi;
    iter=iter+1;
    t=t+dt;
    s=s+1;
    
    %vorticity
    for i=1:n
        for j=1:n
            if j~=1 && j~=n && i~=1 && i~=n
                
                uu=(psi(i,j+1)-psi(i,j-1));
                umin=min(uu,0);uplus=max(uu,0);
                vv=-(psi(i+1,j)-psi(i-1,j));
                vmin=min(vv,0);vplus=max(vv,0);
                
                % 1st-O upwind for blending
                vortxmin= (vort(i,j)-vort(i-1,j))/h;
                vortxplus = (vort(i+1,j)-vort(i,j))/h;
                vortymin= (vort(i,j)-vort(i,j-1))/h;
                vortyplus = (vort(i,j+1)-vort(i,j))/h;
                
                
                
                vort(i,j)=vort(i,j)+dt*...
                    (-uu*(vort(i+1,j)-vort(i-1,j))/(4*h^2)...
                    +q*(uplus*vortxmin + umin*vortxplus)...
                    -vv*(vort(i,j+1)-vort(i,j-1))/(4*h^2)...
                    +q*(vplus*vortymin + vmin*vortyplus)...
                    +nu*(vort(i+1,j)+vort(i,j+1)+vort(i-1,j)+vort(i,j-1)-4*vort(i,j))/h^2);
            elseif j==1 && i~=1 && i~=n
                vort(i,j)=2*(psi(i,j)-psi(i,j+1))/h^2; %Bottom
            elseif j==n && i~=1 && i~=n
                vort(i,j)=2*(psi(i,j)-psi(i,j-1))/h^2-2*U/h; % Top
            elseif i==1 && j~=1 && j~=n
                vort(i,j)=2*(psi(i,j)-psi(i+1,j))/h^2; % Left
            elseif i==n && j~=1 && j~=n
                vort(i,j)=2*(psi(i,j)-psi(i-1,j))/h^2; % Right
            elseif (i==1 && j==1) % Set corners equal to neighbors or zero
                vort(i,j)=0;
            elseif (i==n && j==1)
                vort(i,j)=0;
            elseif (i==1 && j==n)
                vort(i,j)=(vort(i+1,j));
            elseif (i==n && j==n)
                vort(i,j)=(vort(i-1,j));
            end
        end
    end
 
    % streamfunction 
    for i=1:n
        for j=1:n
            if j~=1 && j~=n && i~=1 && i~=n
                psi(i,j)=0.25*(h^2*vort(i,j)+psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1));
            elseif y(j)==0.3 % top
                psi(i,j)=0.07;
            elseif x(i)==0 && y(j)>0.25 && y(j)<=0.3 % left 
                psi(i,j)=0.07;
            elseif x(i)==0.3 && y(j)>.1 && y(j)<=0.3 % right
                psi(i,j)=0.07;
            else
                psi(i,j)=0;
            end
        end
    end
    % calculate velocities
    for i=1:n
        for j=1:n
            if j~=1 && j~=n && i~=1 && i~=n
                u(i,j)=(psi(i,j+1)-psi(i,j-1))/(2*h);
                v(i,j)=(psi(i+1,j)-psi(i-1,j))/(2*h);
            elseif j==n
                u(i,j)=1;
                v(i,j)=0;
            else
                u(i,j)=0;
                v(i,j)=0;
            end
        end
    end
    
%     figure(4);
%     subplot(1,3,1);
%     velocityfig=quiver(x,y,u',v');axis square;axis([0 L 0 L]);
% %     velmagfig=contourf(x,y,sqrt(u'.^2+v'.^2));axis square;axis([0 L 0 L]);
%     subplot(1,3,2);
%     streamfig=contourf(x,y,psi',20,'LineColor','none' );axis square;axis([0 L 0 L]);
%     subplot(1,3,3);
%     vortfig=contourf(x,y,abs(vort'),40,'LineColor','none' );axis square;axis([0 L 0 L]);
%     drawnow;
    
end
elapsedTime=toc;
fprintf('Took %0.4f seconds \n',elapsedTime);
%% plots
figure(1)
quiver(x,y,u',v')
axis([0 L 0 L]);
xlabel('x');ylabel('y');
title(sprintf('Velocity Field Re=200 q=%0.2f t=%0.2f',q,t));

figure(2)
contourf(x,y,psi',20,'LineColor','none' )
axis([0 L 0 L])
xlabel('x');ylabel('y');
colorbar
title(sprintf('Streamlines Re=200 q=%0.2f t=%0.2f',q,t));

figure(3)
contourf(x,y,abs(vort'),40,'LineColor','none' )
axis([0 L 0 L])
xlabel('x');ylabel('y');
colorbar
title(sprintf('Vorticity Re=200 q=%0.2f t=%0.2f',q,t));

%%
% % get values at desired pts
% for kk=1:4
%     inds=[51 251 751 1501];
%     uuuu=uout{inds(kk)};
%     vvvv=vout{inds(kk)};
%     vortt=vortout{inds(kk)};
%     times=[0.1 0.5 1.5 3];
%     t= times(kk);
%     fprintf('x-Velocity (u) @ (0.15,0.25) q=%0.2f t=%0.1f is %0.3f \n',q,t,uuuu(16,26));
%     fprintf('y-Velocity (v) @ (0.15,0.25) q=%0.2f t=%0.1f is %0.3f \n',q,t,vvvv(16,26));
%     fprintf('Vorticity @ (0.15,0.25) q=%0.2f t=%0.1f is %0.3f \n',q,t,vort(26,16));
%     
%     
%     fprintf('x-Velocity (u) @ (0.0,0.22) q=%0.2f t=%0.1f is %0.3f \n',q,t,uuuu(1,23));
%     fprintf('y-Velocity (v) @ (0.0,0.22) q=%0.2f t=%0.1f is %0.3f \n',q,t,vvvv(16,23));
% end
