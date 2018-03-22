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
title('OBJECT + 2*EOBJ','FontSize',14,'FontWeight','Bold');
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
F = ECLK | EOBJ;

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

figure('Color','w');
a = imagesc(xa,ya,F');
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
title('FORCE MAP','FontSize',14,'FontWeight','Bold');
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM #4 -- CALCULATE COORDINATES BY SOLVING LAPLACE’S EQUATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE DERIVATICE MATRICES
NS = [Nx Ny];
RES = [dx dy];
BC = [1 1];
[DX,D2X,DY,D2Y] = fdder(NS,RES,BC);

% BUILD L MATRIX
F = diag(sparse(F(:)));
I = speye(Nx*Ny);
L = D2X + D2Y;
L = F + (I - F)*L;

% REDUCE PROBLEM
M = ECLK | (CLOAK & ~OBJECT) | EOBJ;
ind = find(M(:));
L = L(ind,ind);
XF = XF(ind);
YF = YF(ind);

% SOLVE LAPLACE'S EQUATIONS FOR X-VALUES
u = L\XF;

% REINSERT INTO GRID
XB = zeros(Nx,Ny);
XB(ind) = u;

% SOLVE LAPLACE'S EQUATIONS FOR Y-VALUES
u = L\YF;

% REINSERT INTO GRID
YB = zeros(Nx,Ny);
YB(ind) = u;

% VISUALIZE
figure('Color','w');
subplot(121);
a = imagesc(xa,ya,XB');
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
title('XB','FontSize',14,'FontWeight','Bold');
colorbar;

subplot(122);
a = imagesc(xa,ya,YB');
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
title('YB','FontSize',14,'FontWeight','Bold');
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM #5 -- CALCULATE PERMITTIVITY TENSOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE TENSORS
ERxx = ones(Nx,Ny);
ERxy = zeros(Nx,Ny);
ERxz = ERxy;
ERyx = ERxy;
ERyy = ERxx;
ERyz = ERxy;
ERzx = ERxy;
ERzy = ERxy;
ERzz = ERxx;

% CALCULATE DERIVATIVES OF TRANSFORMED GRID
Dxx = reshape(DX*XB(:),Nx,Ny);
Dxy = reshape(DY*XB(:),Nx,Ny);
Dyx = reshape(DX*YB(:),Nx,Ny);
Dyy = reshape(DY*YB(:),Nx,Ny);
M = (CLOAK & ~OBJECT);

% BUILD ER
for ny = 1:Ny
    for nx = 1:Nx
        if M(nx,ny)
            % BUILD BACKGROUND TENSOR
            ER = [ERxx(nx,ny) ERxy(nx,ny) ERxz(nx,ny);
                ERyx(nx,ny) ERyy(nx,ny) ERyz(nx,ny);
                ERzx(nx,ny) ERzy(nx,ny) ERzz(nx,ny)];
            
            % BUILD JACOBIAN MATRIX
            J = [Dxx(nx,ny) Dxy(nx,ny) 0;
                Dyx(nx,ny) Dyy(nx,ny) 0;
                0          0      1];
            
            % TRANSFORM ER
            J = inv(J);
            ER = J*ER*J.'/det(J);
            
            % POPULATE TENSORS
            ERxx(nx,ny) = ER(1,1);
            ERxy(nx,ny) = ER(1,2);
            ERxz(nx,ny) = ER(1,3);
            ERyx(nx,ny) = ER(2,1);
            ERyy(nx,ny) = ER(2,2);
            ERyz(nx,ny) = ER(2,3);
            ERzx(nx,ny) = ER(3,1);
            ERzy(nx,ny) = ER(3,2);
            ERzz(nx,ny) = ER(3,3);
        end
    end
end

% VISUALIZE
figure('Color','w');
subplot(331);
imagesc(xa,ya,ERxx');
colorbar;
caxis([-10 10]);
axis equal tight;
axis off;
title('ERxx','FontSize',12);
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);

subplot(332);
imagesc(xa,ya,ERxy');
colorbar;
caxis([-10 10]);
axis equal tight;
axis off;
title('ERxy','FontSize',12);
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);

subplot(333);
imagesc(xa,ya,ERxz');
colorbar;
caxis([-10 10]);
axis equal tight;
axis off;
title('ERxz','FontSize',12);
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);

subplot(334);
imagesc(xa,ya,ERyx');
colorbar;
caxis([-10 10]);
axis equal tight;
axis off;
title('ERyx','FontSize',12);
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);

subplot(335);
imagesc(xa,ya,ERyy');
colorbar;
caxis([-10 10]);
axis equal tight;
axis off;
title('ERyy','FontSize',12);
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);

subplot(336);
imagesc(xa,ya,ERyz');
colorbar;
caxis([-10 10]);
axis equal tight;
axis off;
title('ERyz','FontSize',12);
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);

subplot(337);
imagesc(xa,ya,ERzx');
colorbar;
caxis([-10 10]);
axis equal tight;
axis off;
title('ERzx','FontSize',12);
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);

subplot(338);
imagesc(xa,ya,ERzy');
colorbar;
caxis([-10 10]);
axis equal tight;
axis off;
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);

subplot(339);
imagesc(xa,ya,ERzz');
colormap('hot');
colorbar;
caxis([-5 5]);
axis equal tight;
axis off;
title('ERzz','FontSize',12);
xlabel('x','FontSize',12);
ylabel('y','FontSize',12,'Rotation',0);