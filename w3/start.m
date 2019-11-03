close all;
clearvars;
clc

%I=double(imread('zigzag_mask.png'));
%I=mean(I,3); %To 2D matrix
%I=double(imread('circles.png'));
%I=double(imread('noisedCircles.tif'));
%I=double(imread('phantom17.bmp'));
%I=double(imread('phantom18.bmp'));
%I=double(imread('wrench.bmp'));
%I=double(imread('camera.png'));
%I=double(imread('cell.png'));
I=double(imread('coins.png'));

%Normalization of the image to [0 1]
I=mean(I,3);
I=I-min(I(:));
I=I/max(I(:));

[ni, nj]=size(I);

%Lenght and area parameters
%circles.png mu=1, mu=2, mu=10
%noisedCircles.tif mu=0.1
%phantom17 mu=1, mu=2, mu=10
%phantom18 mu=0.2 mu=0.5
%hola carola mu=1
mu=0.25;
nu=0;

%Averages parameters
%hola carola lambda1=10^-3, lambda2=10^-3
lambda1=1;
lambda2=1;

epHeaviside=1;
eta=1e-8;
tol=1e-3;
dt=0.5;
iterMax=10000;
reIni=0;

%%Initial phi
[X, Y]=meshgrid(1:nj, 1:ni);
%phi_0=(-sqrt( ( X-round(ni/2)).^2 + (Y-round(nj/2)).^2)+50);
phi_0 = sin(pi/5 * X) .* sin(pi/5 * Y);
%This initialization allows a faster convergence for phantom 18
%phi_0=(-sqrt( ( X-round(ni/2)).^2 + (Y-round(nj/4)).^2)+50);
%For the Hola carola problem
%phi_0=I;

%Normalization of the initial phi to [-1 1]
phi_0=phi_0-min(phi_0(:));
phi_0=2*phi_0/max(phi_0(:));
phi_0=phi_0-1;

%%Explicit Gradient Descent
seg=sol_ChanVeseIpol_GDExp(I, phi_0, mu, nu, eta, lambda1, lambda2, tol, epHeaviside, dt, iterMax, reIni);

figure;
subplot(1,2,1);
imshow(I);
subplot(1,2,2)
imshow(seg);
