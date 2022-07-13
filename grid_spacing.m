%%ME303 Project 2 Part 1 - Jacob Chartrand, Evan Bernat, Jared Elliott,
%%Chris Gosk
clear all
close all
clc
format long
%Assumptions
k = 1; %Thermal conductivity
L = 1;
M = 30000; %Time resolution
time_res_lim = 45;
c = hot(time_res_lim);

Temp_diff = [];
old_Temp = zeros(1,M);

for w = 1:1:time_res_lim
%clearvars -except w old_Temp Temp_diff c time_res_lim
clear i j T0 T1 Temp
%Initilization
k = 1; %Thermal conductivity
L = 1;
T = 10; %t in (0,T)
N = w; %Space resolution
M = 30000; %Time resolution
dx = L/N; dt = T/M; %Grid spacing
alpha = k*dt/dx^2;
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

p = plot(0:dx:L,Temp(round(0.01*(1/dt)),:));
p.Color = c(time_res_lim - w + 1,:);
set(gca,'Color',[0.7 0.7 0.7])
hold on
drawnow;

Temp_diff = [Temp_diff; mean(old_Temp) - mean(Temp(round(0.01*(1/dt)),:))];

old_Temp = mean(Temp(round(0.01*(1/dt)),:));
end
%% Plotting
figure(2)
semilogy((2:1:N),abs(Temp_diff(2:N)))

xlabel('N value')
ylabel('Difference in Average Temp vs N-1')
xlim([0 time_res_lim + 1])
