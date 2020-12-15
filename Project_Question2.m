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

%% ADI solution

u_ADI = u;
v_ADI = v;
u0_ADI = u;   % for n+1/2
v0_ADI = v;   % for n+1/2

for i=2:N-1
    
    for j=2:N-1
        Aa1 = u_ADI(i,j) * (u0_ADI(i+1,j) - u0_ADI(i-1,j))/(2*dx) + ...
            v_ADI(i,j) * (u_ADI(i,j+1) - u_ADI(i,j-1))/(2*dy);
        Aa2 = nju * (u0_ADI(i+1,j) - 2*u0_ADI(i,j) + u0_ADI(i-1,j))/(dx^2) + ...
            nju * (u_ADI(i,j+1) - 2*u_ADI(i,j) + u_ADI(i,j-1))/(dy^2);       
        Ba1 = u0_ADI(i,j) * (u0_ADI(i+1,j) - u0_ADI(i-1,j))/(2*dx) + ...
            v0_ADI(i,j) * (u_ADI(i,j+1) - u_ADI(i,j-1))/(2*dy);
        Ba2 = nju * (u0_ADI(i+1,j) - 2*u0_ADI(i,j) + u0_ADI(i-1,j))/(dx^2) + ...
            nju * (u_ADI(i+1,j) - 2*u_ADI(i,j) + u_ADI(i,j-1))/(dy^2);                
        Ca1 = u_ADI(i,j) * (v0_ADI(i+1,j) - v0_ADI(i-1,j))/(2*dx) + ...
            v_ADI(i,j) * (v_ADI(i,j+1) - v_ADI(i,j-1))/(2*dy);
        Ca2 = nju * (v0_ADI(i+1,j) - 2*v0_ADI(i,j) + v0_ADI(i-1,j))/(dx^2) + ...
            nju * (v_ADI(i,j+1) - 2*v_ADI(i,j) + v_ADI(i,j-1))/(dy^2);                 
        Da1 = u0_ADI(i,j) * (v0_ADI(i+1,j) - v0_ADI(i-1,j))/(2*dx) + ...
            v0_ADI(i,j) * (v_ADI(i,j+1) - v_ADI(i,j-1))/(2*dy);
        Da2 = nju * (v0_ADI(i+1,j) - 2*v0_ADI(i,j) + v0_ADI(i-1,j))/(dx^2) + ...
            nju * (v_ADI(i,j+1) - 2*v_ADI(i,j) + v_ADI(i,j-1))/(dy^2);           
        
        u_ADI(i,j) = u_ADI(i,j) + dt/2 * (Aa2 - Aa1) + dt/2 * (Ba2 - Ba1);
        v_ADI(i,j) = v_ADI(i,j) + dt/2 * (Ca2 - Ca1) + dt/2 * (Da2 - Da1);
    end
    
end

%% Plotting surface of u(x,y) and v(x,y)
figure(1)
[X,Y] = meshgrid(x,y);
surf(X,Y,u_ADI);hold on;
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
surf(X,Y,v_ADI);hold on;
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