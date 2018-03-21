% Homework
%
% This MATLAB code designs an invisibility cloak using
% transformation electromagnetics.
%
% EE5322 21st Century Electromagnetics
% Spring 2016
%
% Instructor: Dr. Raymond C. Rumpf

% INITIALIZE MATLAB
close all;
clc;
clear all;

% OPEN FIGURE WINDOW
fig = figure('Color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OBJECTS
w = 0.9; %dimension of square cloak
t = 0.4; %length of triangle side

% GRID
Sx = 1;
Sy = Sx;
Nx = 100;
Ny = round(Nx*Sy/Sx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM #1 -- BUILD CLOAK AND OBJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE GRID RESOLUTION
dx = Sx/Nx;
xa = [0 : Nx-1]*dx; xa = xa - mean(xa);
dy = Sy/Ny;
ya = [0 : Ny-1]*dy; ya = ya - mean(ya);
[Y,X] = meshgrid(ya,xa);

% COMPUTE START AND STOP INDICES
x1 = ceil((Sx-w)/2/dx) + 1;
x2 = Nx - x1 + 1;
y1 = ceil((Sy-w)/2/dy) + 1;
y2 = Ny - y1 + 1;

% FILL CLOAK
CLOAK = zeros(Nx,Ny);
CLOAK(x1:x2,y1:y2) = 1;

% VISUALIZE CLOAK
a = imagesc(xa,ya,CLOAK');
a = get(a,'Parent');
set(a,'FontSize',10);
T = [-0.5 : 0.25 : +0.5];
L = {};
for n = 1 : length(T)
    L{n} = num2str(T(n),'%1.2f');
end
L{3} = 0;
set(a,'XTick',T,'XTickLabel',L,'YTick',T,'YTickLabel',L);
colormap('gray');
axis equal tight
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);
title('CLOAK','FontSize',14,'FontWeight','Bold');

% COMPUTE START AND STOP INDICES
OBJECT = zeros(Nx,Ny);
h   = sqrt(t^2 - (t/2)^2);    % Height of triangle
ny  = round(h/dy);
ny1 = 1 + floor((Ny - ny)/2);
ny2 = ny1 + ny - 1;
for ny = ny1 : ny2
    f = (ny - ny1 + 1)/(ny2 - ny1 + 1);
    nx = round(f*t/dx);
    nx1 = 1 + floor((Nx - nx)/2);
    nx2 = nx1 + nx - 1;
    OBJECT(nx1:nx2,ny) = 1;
end 

% VISUALIZE OBJECT
fig2 = figure('Color','w');
a = imagesc(xa,ya,OBJECT');
a = get(a,'Parent');
set(a,'FontSize',10);
T = [-0.5 : 0.25 : +0.5];
L = {};
for n = 1 : length(T)
    L{n} = num2str(T(n),'%1.2f');
end
L{3} = 0;
set(a,'XTick',T,'XTickLabel',L,'YTick',T,'YTickLabel',L);
colormap('gray');
axis equal tight
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);
title('OBJECT','FontSize',14,'FontWeight','Bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM #2 – DETECT EDGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE EDGE ARRAYS
ECLK = zeros(Nx,Ny);
EOBJ = ECLK;

% DETECT EDGES
for ny = 2:Ny-1
    for nx = 2:Nx-1
        if ~CLOAK(nx,ny) & (sum(sum(CLOAK(nx-1:nx+1,ny-1:ny+1))) > 0)
            ECLK(nx,ny) = 1;
        end
        if OBJECT(nx,ny) & (sum(sum(OBJECT(nx-1:nx+1,ny-1:ny+1))) < 8)
            EOBJ(nx,ny) = 1;
        end
    end
end

% VISUALIZE
figure('Color','w');
a = imagesc(xa,ya,(2*EOBJ + OBJECT)');
a = get(a,'Parent');
set(a,'FontSize',10);
T = [-0.5 : 0.25 : +0.5];
L = {};
for n = 1 : length(T)
    L{n} = num2str(T(n),'%1.2f');
end
L{3} = 0;
set(a,'XTick',T,'XTickLabel',L,'YTick',T,'YTickLabel',L);
colormap('gray');
axis equal tight
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);
title('OBJECT + 2*ECLK','FontSize',14,'FontWeight','Bold');
colorbar;

figure('Color','w');
a = imagesc(xa,ya,(2*ECLK + CLOAK)');
a = get(a,'Parent');
set(a,'FontSize',10);
T = [-0.5 : 0.25 : +0.5];
L = {};
for n = 1 : length(T)
    L{n} = num2str(T(n),'%1.2f');
end
L{3} = 0;
set(a,'XTick',T,'XTickLabel',L,'YTick',T,'YTickLabel',L);
colormap('gray');
axis equal tight
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);
title('CLOAK + 2*ECLK','FontSize',14,'FontWeight','Bold');
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM #3 -- GENERATE BOUNDARY CONDITIONS FOR SPATIAL TRANSFORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FORCE MAP
F = ECLK + EOBJ;

% GENERATE INITIAL BOUNDARY VALUES
XA = X.*F;
YA = Y.*F;

% GENERATE FORCED COORDINATES
XF = X.*ECLK;
YF = Y.*ECLK;

% VISUALIZE
figure('Color','w');
subplot(221);
a = imagesc(xa,ya,XA');
a = get(a,'Parent');
set(a,'FontSize',10);
T = [-0.5 : 0.25 : +0.5];
L = {};
for n = 1 : length(T)
    L{n} = num2str(T(n),'%1.2f');
end
L{3} = 0;
set(a,'XTick',T,'XTickLabel',L,'YTick',T,'YTickLabel',L);
colormap('hot');
axis equal tight
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);
title('XA','FontSize',14,'FontWeight','Bold');
colorbar;

subplot(222);
a = imagesc(xa,ya,YA');
a = get(a,'Parent');
set(a,'FontSize',10);
T = [-0.5 : 0.25 : +0.5];
L = {};
for n = 1 : length(T)
    L{n} = num2str(T(n),'%1.2f');
end
L{3} = 0;
set(a,'XTick',T,'XTickLabel',L,'YTick',T,'YTickLabel',L);
colormap('hot');
axis equal tight
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);
title('YA','FontSize',14,'FontWeight','Bold');
colorbar;

subplot(223);
a = imagesc(xa,ya,XF');
a = get(a,'Parent');
set(a,'FontSize',10);
T = [-0.5 : 0.25 : +0.5];
L = {};
for n = 1 : length(T)
    L{n} = num2str(T(n),'%1.2f');
end
L{3} = 0;
set(a,'XTick',T,'XTickLabel',L,'YTick',T,'YTickLabel',L);
colormap('hot');
axis equal tight
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);
title('XF','FontSize',14,'FontWeight','Bold');
colorbar;

subplot(224);
a = imagesc(xa,ya,YF');
a = get(a,'Parent');
set(a,'FontSize',10);
T = [-0.5 : 0.25 : +0.5];
L = {};
for n = 1 : length(T)
    L{n} = num2str(T(n),'%1.2f');
end
L{3} = 0;
set(a,'XTick',T,'XTickLabel',L,'YTick',T,'YTickLabel',L);
colormap('hot');
axis equal tight
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);
title('YF','FontSize',14,'FontWeight','Bold');
colorbar;
