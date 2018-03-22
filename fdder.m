function [DX,D2X,DY,D2Y] = fdder(NS,RES,BC)
% FDDDER Finite-Difference Derivative Operators
%
% [DX,D2X,DY,D2Y] = fdder(NS,RES,BC);
%
% This MATLAB function generates matrix derivative operators
% for scalar functions on a collocated grid.
%
% INPUT ARGUMENTS
% ================
% NS [Nx Ny] size of grid
% RES [dx dy] grid resolution
% BC [BCx BCy] Boundary Conditions
% 0=Dirichlet, -1=Periodic, +1=Neumann
%
% OUTPUT ARGUMENTS
% ================
% DX First-order derivative with respect to x
% D2X Second-order derivative with respect to x
% DY First-order derivative with respect to y
% D2Y Second-order derivative with respect to y

% Misc. Housekeeping
Nx  = NS(1);
Ny  = NS(2);
dx  = RES(1);
dy  = RES(2);
BCx = BC(1);
BCy = BC(2);

% Size of Matrices
M = Nx*Ny;

% Utility Matrix
Z = sparse(M,M);
O = ones(Nx,Ny);
O = O(:);

% Handle Special Cases
if (Nx == 1)
  DX  = Z;
  D2X = Z;
elseif (Nx ~= 1)
  DX  = Z;
  D2X = Z;
  DX  = spdiags(-1*O,-1,DX);
  DX  = spdiags(+1*O,+1,DX);
  D2X = spdiags(O,-1,D2X);
  D2X = spdiags(-2*O,0,D2X);
  D2X = spdiags(O,+1,D2X);
end 

if (Ny == 1)
  DY  = Z;
  D2Y = Z;
elseif (Ny ~= 1)
  DY  = Z;
  D2Y = Z;
  DY  = spdiags(+1*O,Nx,DY);
  DY  = spdiags(-1*O,-Nx,DY);
  D2Y = spdiags(O,+Nx,D2Y);
  D2Y = spdiags(O,-Nx,D2Y);
  D2Y = spdiags(-2*O,0,D2Y);
end

%
%% IMPLEMENT BOUNDARY CONDITIONS
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dirichlet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DIRICHLET ALONG X (DY and D2Y stay the same for Dirichlet)
if(BCx == 0) 
  for ny = 1 : Ny-1 
    m = Nx + (ny-1)*Nx;
    n = m + 1;
    DX(m,n)   = 0;
    DX(n,m)   = 0;
    D2X(m,n)  = 0;
    D2X(n,m)  = 0;
  end 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERIODIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PERIODIC ALONG X
if (BCx == -1 && Nx ~= 1)
  for ny = 1 : Ny
    m = Nx + (ny-1)*Nx;
    n = m - (Nx-1);
    DX(m,n)   = 1;
    DX(n,m)   = -1;
    D2X(m,n)  = 1;
    D2X(n,m)  = 1;
  end 
  for ny = 1 : Ny-1
    m = Nx + (ny-1)*Nx;
    n = m + 1;
    DX(m,n)   = 0;
    DX(n,m)   = 0;
    D2X(m,n)  = 0;
    D2X(n,m)  = 0;
  end 
end 

% PERIODIC ALONG Y
if(BCy == -1 && Ny ~= 1)
  for nx = 1 : Nx 
    m = Nx*Ny - (nx-1);
    n = Nx - (nx-1);
    DY(m,n)   = 1;
    DY(n,m)   = -1;
    D2Y(m,n)  = 1;
    D2Y(n,m)  = 1;
  end 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NEUMANN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NEUMANN ALONG X
if(BCx == +1 && Nx ~= 1) 
  m = 1;
  n = 1;
  for ny = 1 : Ny 
    DX(m,n)     = -2;
    DX(m,n+1)   = +2;
    D2X(m,n)    = 0;
    D2X(m,n+1)  = 0;
    m           = Nx*ny;
    n           = Nx*ny-1;
    DX(m,n)     = -2;
    DX(m,n+1)   = +2;
    D2X(m,n)    = 0;
    D2X(m,n+1)  = 0;
    m           = m + 1;
    n           = n + 2;
  end 
  
  for ny = 1 : Ny-1
    m         = Nx + (ny-1)*Nx;
    n         = m + 1;
    DX(m,n)   = 0;
    DX(n,m)   = 0;
    D2X(m,n)  = 0;
    D2X(n,m)  = 0;
  end 
end 

% NEUMANN ALONG Y
if(BCy == +1 && Ny ~= 1)
  for nx = 1 : Nx
    m = nx;
    n = Nx*Ny-nx+1;
    o = m+Nx;
    p = n-Nx;
    DY(m,m)   = -2;
    DY(n,n)   = +2;
    DY(m,o)   = +2;
    DY(n,p)   = -2;
    D2Y(m,m)  = 0;
    D2Y(n,n)  = 0;
    D2Y(m,o)  = 0;
    D2Y(n,p)  = 0;
  end 
end 

% FINISH DERIVATIVE MATRICES
DX  = DX/(2*dx);
D2X = D2X/(dx^2);
DY  = DY/(2*dy);
D2Y = D2Y/(dy^2);

end