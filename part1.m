%%ME303 Project 2 Part 1 - Jacob Chartrand, Evan Bernat, Jared Elliott,
%%Chris Gosk
clear all
close all
clc
k = 0.00991; %thermal conductivity of egg white
%k = 1;
R = 1; %x in (0,R)
T = 50; %t in (0,T)
N = 20; %Space resolution
M = 5000; %Time resolution
dx = R/N; dt = T/M; %Grid spacing
alpha = k*dt/dx^2;
temp_w = 100; %Water temperature
temp_egg = 25; %Inital egg temperature

%Node Position
for i = 1:N+1
x(i) = (i-1)*dx;
end

%IC
for i = 1:N+1
%T0(i) = sin(pi*x(i));
T0(i) = cos(x(i)-pi/2)
end

%Explicit PDE 
for j = 1:M
for i = 2:N
T1(i) = T0(i) + alpha*(T0(i+1)-2*T0(i)+T0(i-1));
end

T1(1) = 0;
T1(N+1) = 0;

T0 = T1;
Temp(j,:) = T1;
end

%% Plotting
[X,Y] = meshgrid(0:dx:R,dt:dt:T);
mesh(X,Y,Temp); colormap('hot');
xlabel('x')
ylabel('t')
zlabel('T(x,t')
colorbar

