% Gauss Seidel w/SOR to solve Laplace Eqn
%     T_xx + T_yy = 0
% with Neumann and Dirichlet BCs.
%%
clearvars;close all;
% Parameters
nx=21; ny=21;                  % number of space steps
x=linspace(0,2,nx);            % x range
y=linspace(0,1,ny);            % y range
dx=x(2)-x(1);dy=y(2)-y(1);
B=(dx/dy)^2;
T_gs=zeros(nx,ny);             % initialize T for Point Jacobi
tol=1e-5;                      % error tolerance
err=1;                         % error
k=1;                           % iteration index
omega=1.95;                    % (optimal) relaxation parameter
%%
% Boundary Conditions
T_gs(1,:)=0;%left
T_gs(nx,:)=y;%right
T_gs(:,1)=T_gs(:,2);%bottom
T_gs(:,end)=T_gs(:,end-1);%top

Texact=T_gs;     %analytical 
%%
% Time loop
i=1;j=1;
rms=0;
while err>tol
    T_gsold=T_gs;

    for i=2:nx-1
        for j=2:ny-1
            %Gauss-Seidel
            T_gs(i,j)= (1-omega)*T_gsold(i,j) +(omega/(2*(1+B)))*(T_gs(i-1,j)+T_gsold(i+1,j)+B*(T_gs(i,j-1)+T_gsold(i,j+1)));
        end
    end

    %boundary conditions
    T_gs(1,:)=0;%left
    T_gs(nx,:)=y;%right
    T_gs(:,1)=T_gs(:,2);%bottom
    T_gs(:,end)=T_gs(:,end-1);%top
    
    err= max(max(abs(T_gs-T_gsold)));
    k=k+1;

end
%%
% analytical solution
for iii=1:nx
    for jjj=1:ny  
        A=0;
        for n=1:101
            if mod(n,2)==1
                A = A +((n*pi)^-2 * csch(2*n*pi) * sinh(n*pi*x(iii)) * cos(n*pi*y(jjj)));
            end 
        end
        Texact(iii,jjj)=(x(iii)/4)-(4*A);
    end
end
%%
levs=linspace(0,.95,16);
contour(x,y,T_gs',levs);
hold on
contour(x,y,Texact',levs,'--');
xlabel('x');ylabel('y');
axis square;
title('T(x,y) Gauss-Seidel (w/ SOR)');
legend('Gauss-Seidel','Analytical','location','best');
colorbar

pts={[0.5 0.25] [0.5 0.75] [1.5 0.25] [1.5 0.75]};
for zz=1:4
    pt=pts{zz};
    ii=find(x==pt(1));jj=find(y==pt(2));
    fprintf('T(%.1f,%.2f)= %f',pt(1),pt(2),T_gs(ii,jj)); fprintf('\n');
    fprintf('Texact(%.1f,%.2f)= %f',pt(1),pt(2),Texact(ii,jj)); fprintf('\n');
end   
fprintf('%i iterations \n',k);

%rel root mean square error 
tf = false(size(T_gs));
tf(2:end-1,2:end-1) = true ;
rrmse= sqrt(mean((T_gs(tf)./Texact(tf) -1).^2));


