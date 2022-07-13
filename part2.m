%%ME303 Project 2 Part 1 - Jacob Chartrand, Evan Bernat, Jared Elliott,
%%Chris Gosk
clear all
close all
clc

%Assumptions
L = 1;

%Initilization
T = 10; %t in (0,T)
N = 20; %Space resolution
M = 10000; %Time resolution
dx = L/N; dt = T/M; %Grid spacing
alpha = dt/dx^2;

%Node Position
for i = 1:N+1
x(i) = (i-1)*dx;
end

%IC
for i = 1:N+1                                       
T0(i) = cos(pi*x(i));
end

%Explicit method PDE solving 
for j = 1:M %Time
for i = 2:N %Space
T1(i) = T0(i) + alpha*(T0(i+1)-2*T0(i)+T0(i-1));
end

T1(1) = 0;
T1(N+1) = 2;

T0 = T1;
Temp(j,:) = T1;
end

%% Plotting
figure(1)
plot(0:dx:L,Temp(round(0.001*(1/dt)),:))
hold on
plot(0:dx:L,Temp(round(0.01*(1/dt)),:))
hold on
plot(0:dx:L,Temp(round(0.1*(1/dt)),:))
hold on
plot(0:dx:L,Temp(round(10*(1/dt)),:))
hold on
legend('t=0.001','t=0.01','t=0.1','t=10','Location','northwest')

figure(2)
for i = 1:1:M
plot(0:dx:L,Temp(i,:))
title('Temperature vs Location at Varying Times')
xlabel('L (unitless)')
ylabel('Temperature (unitless)')
xlim([0 1]) 
ylim([-2 2])

drawnow;
end
