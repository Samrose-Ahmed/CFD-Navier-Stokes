% Solves 2D Heat eqn
%     T_t = alpha*(T_xx +T_yy)
% with Crank-Nicolson scheme.
%%
clearvars;close all;clc;
tic
%%  parameters
dx = 0.1;               % x space step
x = 0:dx:3.5;           % x range
y = x;                  % y range 
dy = dx;                % y space step
nx = length(x);         % number of space steps
n= (nx-2)^2;            % number of interior nodes
T = zeros(n,1);         % intialize T

alpha=.654;             % diffusion constant
dt = .001;              % time step  
t=0;                    % initialize t
tmax=0.5;               % final t
B = (alpha*dt)/(2*dx^2);
BC=zeros(n,1);         % initialize BC vector
%% BCs

B1=makeMat(nx-2,-B,1-4*B,-B);
A1=makeMat(nx-2,0,-B,0);C1=A1;
Q=makeblkdiag(nx-2,A1,B1,C1);

B2=makeMat(nx-2,B,1+4*B,B);
A2=makeMat(nx-2,0,B,0);C2=A2;
RHS=makeblkdiag(nx-2,A2,B2,C2);

%BCs
Tleft=200;
Tbot=200;
Ttop=0;
Tright=0;
for i2=1:nx-2
    left=(1+(i2-1)*(nx-2));
    right=left+(nx-3);
    BCsec= BC(left:right);
    if(i2==1)
        BCsec(1)=B*(Tbot+Tleft);
        BCsec(2:end-1)=B*Tbot;
        BCsec(end)=B*(Tbot+Tright);
    elseif(i2<nx-2)
        BCsec(1)=B*(Tleft);
        BCsec(2:end-1)=0;
        BCsec(end)=B*(Tright);
    else
        BCsec(1)=B*(Ttop+Tleft);
        BCsec(2:end-1)=B*Ttop;
        BCsec(end)=B*(Ttop+Tright);
    end
    BC(left:right)=BCsec;
end
%% time loop
while t<tmax
    disp(t)
    
    T= Q\((RHS*T)+(2*BC));
    
    t=t+dt;
end

elapsedTime=toc;
disp(elapsedTime);
%% plot
T=reshape(T,[nx-2,nx-2])';
Tmat=zeros(nx);
Tmat(2:end-1,2:end-1)=T;
Tmat(:,1)=Tleft;
Tmat(1,:)=Tbot;

%% functions
function [Output] = makeblkdiag(n,A,B,C)
% function that makes a TDM for Q
% combines multiple matrix A,B,C into one TDM
    Output = zeros(n^2,n^2);
    for i = 1:n:n^2
        for j = i:n:n^2
            Output(i:i+n-1,j:j+n-1) = B;
            break;
        end
    end
    for i = 1:n:n^2-n
        for j = i:n:n^2-n
            x = i+n;
            Output(x:x+n-1,j:j+n-1) = A;
            break;
        end
    end
    for i = 1:n:n^2-n
        for j = i:n:n^2-n
            y = i+n;
            Output(i:i+n-1,y:y+n-1) = C;
            break;
        end
    end
end
function [Output] = makeMat(n,a,b,c)
    % This function makes LHS Matrix 
    Output = zeros(n,n);
    for i = 1:(n+1):(n^2)
        Output(i) = b;
    end
    for i = 2:(n+1):(n^2)
        Output(i) = a;
    end
    for i = (n+1):(n+1):(n^2)
        Output(i) = c;
    end
end
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