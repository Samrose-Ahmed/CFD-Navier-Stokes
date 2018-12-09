% FTCS explicit method for 2D Burgers' eqn
%    phi_t + u*phi_x + v*phi_y = mu(phi_xx + phi_yy)
% uses 1st-O upwind for advection terms and 2nd-O CD for diffusion terms
%%
clearvars
%Parameters
dt=.0005;                        %Width of each time step
it=0;                            %time iteration
t=dt*it;                         %time
dx=.01;                          %Width of space step (x)
dy=dx;                           %Width of space step (y)
x=-1:dx:1;                       %Range of x 
y=0:dy:1;                        %Range of y
nx=length(x);                    %Number of steps in space (x)
ny=length(y);                    %Number of steps in space (y)  
mu=1e-2;                         %Diffusion coefficient
r=(mu*dt)/(dx*dx);               %diffusion courant number
a=10;                            %alpha

phi=zeros(nx,ny);                %Preallocate phi
phi_old=zeros(nx,ny);            %Preallocate phi_old

tol=10^-3;                       %tolerance for error
err=1;                           %error

%%
%Boundary conditions
phi(1,:)=1-tanh(a);%left
phi(nx,:)=1-tanh(a);%right
phi(:,ny)=1-tanh(a);%top
phi(1:round(nx/2),1)=1+tanh((2*x(1:round(nx/2))+1)*a);%bottom

%%
%Calculating at each time step
% i=2:nx-1;j=2:ny-1;
ll=0;

while err>tol
    
%     if mod(t,.1)==0
%         h=contour(x,y,phi','fill','on');
%         colorbar
%         axis([-1 1 0 1])
%         axis square
%         tstr=sprintf('%.4f',t);
%         title({strcat('2-D Burgers'' equation with \mu = ',num2str(mu),', dx=dy=',num2str(dx)), strcat('\phi(t,x,y), time(\itt) = ',tstr,'s')});
%         drawnow;
%         refreshdata(h)
%     end

    phi_old=phi;
    
    for i=2:nx-1
        for j=2:ny-1
            
            u=2*y(j)*(1-x(i)^2);%x velocity
            v= -2*x(i)*(1-y(j)^2);%y velocity
            coax=u*dt/dx;%advection courant number in x
            coay=v*dt/dy;%advection courant number in y
            
            if v>0
                phi(i,j)=phi_old(i,j)...
                    -(coax*(phi_old(i,j)-phi_old(i-1,j))) ... advection in x
                    -(coay.*(phi_old(i,j)-phi_old(i,j-1))) ... advection in y
                    +(r*(phi_old(i+1,j)-4*phi_old(i,j)+phi_old(i-1,j)+phi_old(i,j-1)+phi_old(i,j+1))); % diffusion
            else
                phi(i,j)=phi_old(i,j)...
                    -(coax.*(phi_old(i,j)-phi_old(i-1,j))) ... advection in x
                    -(coay.*(phi_old(i,j+1)-phi_old(i,j))) ... advection in y
                    +(r*(phi_old(i+1,j)-4*phi_old(i,j)+phi_old(i-1,j)+phi_old(i,j-1)+phi_old(i,j+1))); % diffusion
            end     
        end
    end


    %Boundary Conditions
    phi(1,:)=1-tanh(a);%left
    phi(nx,:)=1-tanh(a);%right
    phi(:,ny)=1-tanh(a);%top
    phi(1:round(nx/2),1)=1+tanh((2*x(1:round(nx/2))+1)*a);%bottom 
    phi(round(nx/2)+1:end,1)=phi(round(nx/2)+1:end,2);
    
    err= sum(sum(phi-phi_old)); %error 
    
    %display phi(0,0.5)
    if((ismember(t,[0.5 1 1.5 2])))
        ii=find(x==0);jj=find(y==0.5);
        fprintf('phi(0,0.5) at %.1fs is %f', t, phi(ii,jj)); fprintf('\n');
    end
    if(err<tol && ll==0)
        %plot steady state 
        contour(x,y,phi','fill','on');
        colorbar
        axis([-1 1 0 1])
        axis square
        tstr=sprintf('%.4f',t);
        title({strcat('2-D Burgers'' equation with \mu = ',num2str(mu),', dx=dy=',num2str(dx)), strcat('\phi(t,x,y), time(\itt) = ',tstr,'s')});
        
        %display phi(0,0.5) steady state
        ii=find(x==0);jj=find(y==0.5);
        fprintf('phi(0,0.5) at %1.2fs (steady) is %f', t, phi(ii,jj));fprintf('\n');
        ll=1;
    end
     
    %disp(err);
    it=it+1;
    t=dt*it;
end


