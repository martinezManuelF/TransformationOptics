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
dx = min([w t])/Nx;
xa = [0 : Nx-1]*dx; xa = xa - mean(xa);
dy = Sy/Ny;
ya = [0 : Ny-1]*dy; ya = ya - mean(ya);

% COMPUTE START AND STOP INDICES

