% Code to read the output data file from the Simulink recorder
% Imperial College London
% Flight Sim Lab 18/02/2025
% Author Himmat Kaul
%% read data
clear; clc; close all;

format long
s = settings;
s.matlab.fonts.editor.code.Name.TemporaryValue = 'Calibri';
set(groot,'defaultLineLineWidth',2)  %sets graph line width as 2
set(groot,'defaultAxesFontSize',24)  %sets graph axes font size as 18
set(groot,'defaulttextfontsize',24)  %sets graph text font size as 18
set(groot,'defaultLineMarkerSize',14) %sets line marker size as 8
set(groot,'defaultAxesXGrid','on')   %sets X axis grid on 
set(groot,'defaultAxesYGrid','on')   %sets Y axis grid on
set(groot,'DefaultAxesBox', 'on')   %sets Axes boxes on

picturewidth = 30; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
%%


%%%%%INSERT THE FILENAME HERE%%%%%
fname = "E:\Y2\Flight Dynamics\Group 07 Part2 Lab\FlightDynLabPt2_2025-03-03_11-17-43.mat"; % Update to match the provided file
load(fname) % Load the .mat file

% Check if the loaded data structure exists and extract appropriately
if exist('data', 'var')
    % && isfield(data, 'Data') && ~isempty(data.Data) 
    % Extract all the data
    vInd_kias   = data.Data(:,1); % Indicated Airspeed - Knots
    vTrue_ktas  = data.Data(:,2); % True Airspeed      - Knots
    climb_rate  = data.Data(:,3); % Climb rate         - fpm
    q           = data.Data(:,7); % Pitch rate         - rad/s  
    p           = data.Data(:,8); % Roll rate          - rad/s
    r           = data.Data(:,9); % Yaw rate           - rad/s
    pitch       = data.Data(:,10);% Pitch angle        - deg
    roll        = data.Data(:,11);% Roll angle         - deg
    heading_true= data.Data(:,12);% Heading            - deg
    alpha       = data.Data(:,13);% Angle of Attack    - deg
    beta        = data.Data(:,14);% Side slip angle    - deg
    latitude    = data.Data(:,15);% Latitude           - deg
    longitude   = data.Data(:,16);% Longitude          - deg  
    altitude    = data.Data(:,17);% Altitude           - ft 
    x           = data.Data(:,18);
    y           = data.Data(:,19);
    z           = data.Data(:,20);
    throttle_cmd = data.Data(:,21);% Throttle command  - %
    throttle_actual = data.Data(:,22);
    eng_power   = data.Data(:,23); % Engine power      - hp
    w_empty     = data.Data(:,24); 
    w_payld     = data.Data(:,25); % Payload weight    - lb
    w_fuel      = data.Data(:,26); % Fuel weight       - lb 
    time        = data.Data(:,27); % Time from start   - seconds
    elevator    = data.Data(:,28); % Tail incidence    - % 
    aileron     = data.Data(:,29); %                  - %
    rudder      = data.Data(:,30); %                  - %

    % Plot Time vs True Airspeed to verify results present
    % figure;
    % plot(time, vTrue_ktas, 'b', 'LineWidth', 1.5)
    % xlabel('Time (seconds)')
    % ylabel('True Airspeed (knots)')
    % title('Time vs True Airspeed History')
    % grid on;
    % 
    figure;
    plot(time, altitude, 'b', 'LineWidth', 1.5)
    xlabel('Time (seconds)')
    ylabel('Altitude (feet)')
    title('Time vs Altitude History')
    % grid on;
    % 
    % figure;
    % plot(time, pitch, 'b', 'LineWidth', 1.5)
    % xlabel('Time (seconds)')
    % ylabel('pitch (deg)')
    % title('pitch vs Altitude History')
    % grid on;
    % 
    % figure;
    % plot(time, roll, 'b', 'LineWidth', 1.5)
    % xlabel('Time (seconds)')
    % ylabel('roll (deg)')
    % title('pitch vs roll History')
    % grid on;
    % 
    % figure;
    % plot(time, beta, 'b', 'LineWidth', 1.5)
    % xlabel('Time (seconds)')
    % ylabel('yaw (deg)')
    % title('Yaw vs roll History')
    grid on;
else
    error('The loaded file does not contain the expected data structure.')
end

%%

% convert values to SI units

vInd_kias   = vInd_kias * 0.514444; % Indicated Airspeed - m/s
vTrue_ktas  = vTrue_ktas * 0.514444; % True Airspeed      - m/s
climb_rate  = climb_rate * 0.00508; % Climb rate         - fpm
% q           = data.Data(:,7); % Pitch rate         - rad/s  
% p           = data.Data(:,8); % Roll rate          - rad/s
% r           = data.Data(:,9); % Yaw rate           - rad/s
pitch       = deg2rad(pitch);% Pitch angle        - rad
roll        = deg2rad(roll);% Roll angle         - rad
heading_true= deg2rad(heading_true);% Heading            - rad
alpha       = deg2rad(alpha);% Angle of Attack    - rad
beta        = deg2rad(beta);% Side slip angle    - rad
latitude    = deg2rad(latitude);% Latitude           - rad
longitude   = deg2rad(longitude);% Longitude          - rad  
altitude    = altitude * 0.3048;% Altitude           - m 
% x           = data.Data(:,18);
% y           = data.Data(:,19);
% z           = data.Data(:,20);
% throttle_cmd = data.Data(:,21);% Throttle command  - %
% throttle_actual = data.Data(:,22);
eng_power   = eng_power * 745.7; % Engine power      - W
w_empty     = w_empty * 0.453592; 
w_payld     = w_payld * 0.453592; % Payload weight    - kg
w_fuel      = w_fuel * 0.453592; % Fuel weight       - kg 
% time        = data.Data(:,27); % Time from start   - seconds
% elevator    = data.Data(:,28); % Tail incidence    - % 
% aileron     = data.Data(:,29); %                  - %
% rudder      = data.Data(:,30); %                  - %



%% Testing Functions

% Example format for analysis


start = 2422;
end_time = 2561;

[ ~, idx1 ] = min( abs( time-start ) );
[ ~, idx2 ] = min( abs( time-end_time ) );
% define the start and end points for the sureve from the full data series

% defines "x and y" values to compare example roll subsidance will have roll rate vs time 
V_test = altitude(idx1:idx2);
t_test = time(idx1:idx2);

% Test Figure
figure;
% plots data series
plot(t_test, V_test, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks and Troughs
% [~,id] = findpeaks(V_test, 'MinPeakDistance', 5);
% [~,id2] = findpeaks(-V_test, 'MinPeakDistance', 5);
% displays peaks
%plot(t_test(id), V_test(id1), 'rx',MarkerSize=20)
% plot(t_test(id2), V_test(id2), 'ro',MarkerSize=20)

% labels ect boring stuff
xlabel('Time (seconds)')
ylabel('Altitude')
title('Time vs True Airspeed History - test')
grid on;
hold off



%% Testing Functions

% Example format for analysis

[ ~, idx1 ] = min( abs( time-start ) );
[ ~, idx2 ] = min( abs( time-end_time ) );
% define the start and end points for the sureve from the full data series

% defines "x and y" values to compare example roll subsidance will have roll rate vs time 
V_test = climb_rate(idx1:idx2);
t_test = time(idx1:idx2);

% Test Figure
figure;
% plots data series
plot(t_test, V_test, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks and Troughs
% [~,id] = findpeaks(V_test, 'MinPeakDistance', 5);
% [~,id2] = findpeaks(-V_test, 'MinPeakDistance', 5);
% displays peaks
% plot(t_test(id), V_test(id), 'rx',MarkerSize=20)
% plot(t_test(id2), V_test(id2), 'ro',MarkerSize=20)

% labels ect boring stuff
xlabel('Time (seconds)')
ylabel('climb')
title('Time vs True Airspeed History - test')
grid on;
hold off


%% Testing Functions

% Example format for analysis

[ ~, idx1 ] = min( abs( time-start ) );
[ ~, idx2 ] = min( abs( time-end_time ) );
% define the start and end points for the sureve from the full data series

% defines "x and y" values to compare example roll subsidance will have roll rate vs time 
V_test = pitch(idx1:idx2);
t_test = time(idx1:idx2);

% Test Figure
figure;
% plots data series
plot(t_test, V_test, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks and Troughs
% [~,id] = findpeaks(V_test, 'MinPeakDistance', 5);
% [~,id2] = findpeaks(-V_test, 'MinPeakDistance', 5);
% displays peaks
% plot(t_test(id), V_test(id), 'rx',MarkerSize=20)
% plot(t_test(id2), V_test(id2), 'ro',MarkerSize=20)

% labels ect boring stuff
xlabel('Time (seconds)')
ylabel('Pitch')
title('Time vs True Airspeed History - test')
grid on;
hold off

%% Testing Functions

% Example format for analysis

[ ~, idx1 ] = min( abs( time-start ) );
[ ~, idx2 ] = min( abs( time-end_time ) );
% define the start and end points for the sureve from the full data series

% defines "x and y" values to compare example roll subsidance will have roll rate vs time 
V_test = vTrue_ktas(idx1:idx2);
t_test = time(idx1:idx2);

% Test Figure
figure;
% plots data series
plot(t_test, V_test, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks and Troughs
% [~,id] = findpeaks(V_test, 'MinPeakDistance', 5);
% [~,id2] = findpeaks(-V_test, 'MinPeakDistance', 5);
% displays peaks
% plot(t_test(id), V_test(id), 'rx',MarkerSize=20)
% plot(t_test(id2), V_test(id2), 'ro',MarkerSize=20)

% labels ect boring stuff
xlabel('Time (seconds)')
ylabel('True Airspeed [ms^-1]')
title('Time vs True Airspeed History - test')
grid on;
hold off



%%