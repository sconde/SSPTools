clear all; close all; clc

addpath('../');
N = 300;

%testing = 'ERK';
testing = 'DIRK';

y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);

nE = 15;    % number of elements
K  = 2;     % degree of accuracy

%xmesh = NDG.Mesh1D('Domain', [-1 1], 'NumberElements', nE,...
%    'SolutionDegree', K, 'QuadratureType', 'LGL');

xmesh = NDG.Mesh1D('Domain', [-1 1], 'NumberElements', nE,...
    'SolutionDegree', K, 'QuadratureType', 'Legendre');

