% Point Jacobi Method to solve Laplace Eqn
%     T_xx + T_yy = 0
% with Neumann and Dirichlet BCs.
%%
clearvars
% Parameters
nx=21; ny=21;                  % number of space steps
x=linspace(0,2,nx);            % x range
y=linspace(0,1,ny);            % y range
dx=x(2)-x(1);dy=y(2)-y(1);
B=(dx/dy)^2;
T_pj=zeros(nx,ny);             % initialize T 
tol=1e-5;                      % error tolerance for convergence
err=1;                         % to check convergence
k=1;                           % iteration index
%%
% Boundary Conditions
T_pj(1,:)=0;%left
T_pj(nx,:)=y;%right
T_pj(:,1)=T_pj(:,2);%bottom
T_pj(:,end)=T_pj(:,end-1);%top

Texact=T_pj;     % initialize analytical solution
%%
while err>tol
    T_pjold=T_pj;
    
    for i=2:nx-1
        for j=2:ny-1
            %Point Jacobi
            T_pj(i,j)=(1/(2*((1+B))))*(T_pjold(i-1,j)+T_pjold(i+1,j)+B*(T_pjold(i,j-1)+T_pjold(i,j+1)));
        end
    end

    %boundary conditions
    T_pj(1,:)=0;%left
    T_pj(nx,:)=y;%right
    T_pj(:,1)=T_pj(:,2);%bottom
    T_pj(:,end)=T_pj(:,end-1);%top
    
   
    err= norm(T_pj(:)-T_pjold(:),Inf); 
    k=k+1;  
end
rmspj=rms(T_pj(:)-T_pjold(:));

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
% plot
figure(2)
levs=linspace(0,.95,16);
contour(x,y,T_pj',levs);
hold on
contour(x,y,Texact',levs,'--');
xlabel('x');ylabel('y');
axis square;
title('T(x,y) Point Jacobi');
legend('Point Jacobi','Analytical','location','sw');
colorbar

pts={[0.5 0.25] [0.5 0.75] [1.5 0.25] [1.5 0.75]};
for zz=1:4
    pt=pts{zz};
    ii=find(x==pt(1));jj=find(y==pt(2));
    fprintf('T(%.1f,%.2f)= %f',pt(1),pt(2),T_pj(ii,jj)); fprintf('\n');
    fprintf('Texact(%.1f,%.2f)= %f',pt(1),pt(2),Texact(ii,jj)); fprintf('\n');
end 
fprintf('%i iterations \n',k);

%rel root mean square error 
tf = false(size(T_pj));
tf(2:end-1,2:end-1) = true ;
rrmse= sqrt(mean((T_pj(tf)./Texact(tf) -1).^2));

