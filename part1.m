%%ME303 Project 2 Part 1 - Jacob Chartrand, Evan Bernat, Jared Elliott,
%%Chris Gosk
clear all
close all
clc

%Assumptions

%k = 0.009911; %Average thermal conductivity of egg, W/(C^-1*cm^-1)

k = 0.006; %Thermal conductivity of water
volume = 1260; %Volume of desired egg, cm^3
R = nthroot(volume/((4/3)*pi),3);

% D = 5.35; %Egg diameter as average of egg width and height, cm
% R = D/2; %x in (0,R), cm

temp_w = 100; %Water temperature
temp_egg_init = 12; %Inital egg temperature

%Initilization
T = 5000; %t in (0,T)
N = 20; %Space resolution
M = 10000; %Time resolution
dx = R/N; dt = T/M; %Grid spacing
alpha = k*dt/dx^2;

%Node Position
for i = 1:N+1
x(i) = (i-1)*dx;
end

%IC
for i = 1:N+1
T0(i) = temp_egg_init;
end

%Explicit method PDE solving 
for j = 1:M %Time
for i = 2:N %Space
T1(i) = T0(i) + alpha*(T0(i+1)-2*T0(i)+T0(i-1));
end

T1(1) = T1(2);
T1(N+1) = 100;

T0 = T1;
Temp(j,:) = T1;
end

%% Finish Time Calculation


finish_temp = find(Temp(:,1)>80); %Vector containing indices where
                                         %midpoint is above 80C

finish_time = ((finish_temp(1)/M)*T)+10 %Calculates done time

%% Plotting
[X,Y] = meshgrid(0:dx:R,dt:dt:T);
mesh(X,Y,Temp); colormap('hot');
xlim([0 R])
ylim([0 (finish_time+300)])
zlim([0 100])
title('Ostrich Egg Cooking Plot')
xlabel('x (cm)')
ylabel('t (s)')
zlabel('Temperature (C)')
colorbar

done_plane = patch([0,0,R,R],[finish_time,finish_time,finish_time,finish_time],[100,0,0,100],'k',"FaceAlpha",'0.5');
view(3);
legend([done_plane],'Done Time' )


