%%ME303 Project 2 Part 1 - Jacob Chartrand, Evan Bernat, Jared Elliott,
%%Chris Gosk
clear all
close all
clc
format long

%% Video Recording
Video = VideoWriter('Grid_spacing_animation'); %open video file
Video.FrameRate = 7;
open(Video) %Begins video container

%Assumptions
L = 1; %x in (0,L)
T = 10; %t in (0,T)

%Initilization
Temp_diff = [];
stability_factor = [];
M = 30000; %Time resolution
old_Temp = zeros(1,M);
space_res_lim = 45; %Limit on space resolution for loop
heatmap_2d = hot(space_res_lim); %Create 3 wide array for line coloring


for w = 1:1:space_res_lim
clear i j T0 T1 Temp

N = w; %Space resolution
dx = L/N; dt = T/M; %Grid spacing
alpha = dt/dx^2;
stability_factor = [stability_factor; 1 - 2*alpha]; %Calculate if method is stable
                                                    %Must be >0, used for troubleshooting code

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

%Configuring Plot
f = figure(1); %Initialize figure

    p = plot(0:dx:L,Temp((0.01*(1/dt)),:)); %Plot Temp at t=0.01

    f.WindowState = 'maximized'; %Maximize figure, for better res.
    p.Color = heatmap_2d(space_res_lim - w + 1,:); %Set line color
    set(gca,'Color',[0.7 0.7 0.7]) %Set background color to grey

    color_bar = colorbar; %Initialize colorbar
    colormap hot %Set colorbar to hot
    color_bar.Ticks = [0 1]; %Set colorbar with 2 ticks
    color_bar.TickLabels = {'N=45','N=1'}; %Label colorbar ticks

    title('Temperature Distribution at t=0.01, varying space resolution')
    xlabel('L (unitless)')
    ylabel('Temperature (unitless)')
    xlim([0 1]) 
    ylim([-2 2])

    hold on
    drawnow;

Temp_diff = [Temp_diff; mean(old_Temp) - mean(Temp(round(0.01*(1/dt)),:))];
old_Temp = mean(Temp(round(0.01*(1/dt)),:));

pause(0.01) %Pause and grab frame
frame = getframe(gcf); %get frame
writeVideo(Video, frame);
 
end
close(Video) %End video
%% Plotting
figure(2)
semilogy((2:1:space_res_lim),abs(Temp_diff(2:N)))
title('Change in Average Temperature at Different N Values')
xlabel('N value')
ylabel('Difference in Average Temp vs N-1')
xlim([0 space_res_lim + 1])

figure(3)
plot((1:1:space_res_lim),stability_factor)
title('Stability Factors at Different N Values')
xlabel('N value')
ylabel('Stability Factor')
xlim([0 space_res_lim + 1])


