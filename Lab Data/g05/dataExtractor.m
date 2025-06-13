%Code to read the output data file from the Simulink recorder
%Imperial College London
%Flight Sim Lab 20/21

clear
clc
close all

%%%%%INSERT THE FILENAME HERE%%%%%
fname = 'pitch.mat';
load(fname) %Load the file name

%Extract all the data
vInd_kias   = data.Data(:,1); %Indicated Airspeed - Knots
vTrue_ktas  = data.Data(:,2); %true Airspeed      - Knots
climb_rate  = data.Data(:,3); %climb rate         - fpm
q           = data.Data(:,7); %Pitch rate         - rad/s  
p           = data.Data(:,8); %Roll rate          - rad/s
r           = data.Data(:,9); %Yaw rate           - rad/s
pitch       = data.Data(:,10);%Pitch angle        - deg
roll        = data.Data(:,11);%Roll angle         - deg
heading_true= data.Data(:,12);%Heading            - deg
alpha       = data.Data(:,13);%Angle of Attack    - deg
beta        =data.Data(:,14); %Side slip angle    - deg
latitude    = data.Data(:,15);%Latitude           - deg
longitude   = data.Data(:,16);%Longitude          - deg  
altitude    = data.Data(:,17);%Altitude           - ft 
x           = data.Data(:,18);
y           = data.Data(:,19);
z           = data.Data(:,20);
throttle_cmd = data.Data(:,21);%Throttle command  - %
throttle_actual = data.Data(:,22);
eng_power   = data.Data(:,23); %Engine power      - hp
w_empty     = data.Data(:,24); 
w_payld     = data.Data(:,25); %Payload weight    - lb
w_fuel      = data.Data(:,26); %Fuel weight       - lb 
time        = data.Data(:,27); %Time from start   - seconds
elevator    = data.Data(:,28); %Tail incidence    - % 
aileron     = data.Data(:,29); %                  - %
rudder      = data.Data(:,30); %                  - %
