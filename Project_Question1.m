%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Final Project: Solution of Burger's Equations
%%%% Student Name: Yangchao Liao
%%%% Student ID.: 1299252
%%%% Department: Civil & Environmental Eng.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
close all;
clc;

%% Initial and boundary conditions
x = 0:0.01:1;
y = 0:0.01:1;
N = length(x);  % here, length(y) = length(x) = N
dx = x(2) - x(1);
dy = y(2) - y(1);
  
u = zeros(N,N);
v = zeros(N,N);

for i=1:N
    u(1,i) = 0;  % u(0,y) = 0
    u(N,i) = 0;  % u(1,y) = 0
    v(i,N) = 0;  % v(x,1) = 0
    v(i,1) = 1;  % v(x,0) = 1
    u(i,1) = sin( 2 * pi * x(i) );  % u(x,0) = sin(2*pi*x)
    u(i,N) = sin( 2 * pi * x(i) );  % u(x,1) = sin(2*pi*x)
    v(1,i) = 1 - y(i);  % v(0,y) = 1-y
    v(N,i) = 1 - y(i);  % v(1,y) = 1-y   
end

nju = 0.015;
dt = 0.001;

%% Explicit Euler solution

u_Euler = u;
v_Euler = v;

for i=2:N-1
    
    for j=2:N-1
        A = u_Euler(i,j) * ( u_Euler(i+1,j) - u_Euler(i-1,j) )/(2 * dx) + ...
            v_Euler(i,j) * ( u_Euler(i,j+1) - u_Euler(i,j-1) )/(2 * dy);
        B = nju * ( u_Euler(i+1,j) - 2 * u_Euler(i,j) + u_Euler(i-1,j) )/(dx.^2) + ...
            nju * ( u_Euler(i,j+1) - 2 * u_Euler(i,j) + u_Euler(i,j-1) )/(dy.^2);     
        u_Euler(i,j) = u_Euler(i,j) + dt * ( B - A );

        C = u_Euler(i,j) * ( v_Euler(i+1,j) - v_Euler(i-1,j) )/(2 * dx) + ...
            v_Euler(i,j) * ( v_Euler(i,j+1) - v_Euler(i,j-1) )/(2 * dy);
        D = nju * ( v_Euler(i+1,j) - 2 * v_Euler(i,j) + v_Euler(i-1,j) )/(dx.^2) + ...
            nju * ( v_Euler(i,j+1) - 2 * v_Euler(i,j) + v_Euler(i,j-1) )/(dy.^2);     
        v_Euler(i,j) = v_Euler(i,j) + dt * ( D - C );
    
    
    end
    
end

%% Plotting surface of u(x,y) and v(x,y)
figure(1)
[X,Y] = meshgrid(x,y);
surf(X,Y,u_Euler);hold on;
shading faceted
colormap(cool)

xlabel('x','FontName','Arial','FontSize',25)
ylabel('y','FontName','Arial','FontSize',25)
zlabel('u(x,y)','FontName','Arial','FontSize',25)

set(gca,'linewidth',1.5,'FontName','Arial','FontSize',25);
set(gcf,'Color','w','Units','inches','position',[0,0,8,6]);
AxesH = gca;InSet = get(AxesH, 'TightInset');
set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
hold off;

figure(2)
[X,Y] = meshgrid(x,y);
surf(X,Y,v_Euler);hold on;
shading faceted
colormap(cool)

xlabel('x','FontName','Arial','FontSize',25)
ylabel('y','FontName','Arial','FontSize',25)
zlabel('v(x,y)','FontName','Arial','FontSize',25)

set(gca,'linewidth',1.5,'FontName','Arial','FontSize',25);
set(gcf,'Color','w','Units','inches','position',[0,0,8,6]);
AxesH = gca;InSet = get(AxesH, 'TightInset');
set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
hold off;

%% Plotting by line style
% figure(3)
% for i=1:N 
%    for j=1:N
%     X1(1,1:N) = x(1,i);
%    end
%  plot(X1(1,:),u_Euler(:,i));hold on
% end
% xlabel('x','FontName','Arial','FontSize',25)
% ylabel('u(x,y)','FontName','Arial','FontSize',25)
%  
% figure(4)
% for i=1:N 
%    for j=1:N
%     Y1(1,1:N) = y(1,i);
%    end
%  plot(Y1(1,:),u_Euler(i,:));hold on
%  end
% xlabel('y','FontName','Arial','FontSize',25)
% ylabel('u(x,y)','FontName','Arial','FontSize',25)
%    
% figure(5)
% for i=1:N 
%    for j=1:N
%     X1(1,1:N) = x(1,i);
%    end
%  plot(X1(1,:),v_Euler(:,i));hold on
% end
% xlabel('x','FontName','Arial','FontSize',25)
% ylabel('v(x,y)','FontName','Arial','FontSize',25)
%  
% figure(6)
% for i=1:N 
%    for j=1:N
%     Y1(1,1:N) = y(1,i);
%    end
%  plot(Y1(1,:),v_Euler(i,:));hold on
%  end
% xlabel('y','FontName','Arial','FontSize',25)
% ylabel('v(x,y)','FontName','Arial','FontSize',25)