% Code to read the output data file from the Simulink recorder
% Imperial College London
% Flight Sim Lab 18/02/2025
% Author Himmat Kaul
%% Set Default Figure Parameters
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

%% read data

%%%%%INSERT THE FILENAME HERE%%%%%
fname = "E:\Y2\Flight Dynamics\Lab Data\Group 07\FlightDynLabPt1_2025-02-17_10-45-14.mat"; % Update to match the provided file
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
    aileron     = data.Data(:,29); %                   - %
    rudder      = data.Data(:,30); %                   - %

    % Plot Time vs True Airspeed to verify results present
    figure;
    plot(time, vTrue_ktas, 'b', 'LineWidth', 1.5)
    xlabel('Time (seconds)')
    ylabel('True Airspeed (knots)')
    title('Time vs True Airspeed History')
    grid on;

    figure;
    plot(time, altitude, 'b', 'LineWidth', 1.5)
    xlabel('Time (seconds)')
    ylabel('Altitude (feet)')
    title('Time vs Altitude History')
    grid on;

    figure;
    plot(time, pitch, 'b', 'LineWidth', 1.5)
    xlabel('Time (seconds)')
    ylabel('pitch (deg)')
    title('pitch vs Altitude History')
    grid on;

    figure;
    plot(time, roll, 'b', 'LineWidth', 1.5)
    xlabel('Time (seconds)')
    ylabel('roll (deg)')
    title('pitch vs roll History')
    grid on;

    
else
    error('The loaded file does not contain the expected data structure.')
end

%%

% convert values to SI units

vInd_kias   = vInd_kias * 0.514444; % Indicated Airspeed    - m/s
vTrue_ktas  = vTrue_ktas * 0.514444; % True Airspeed        - m/s
climb_rate  = climb_rate * 0.00508; % Climb rate            - fpm
% q           = data.Data(:,7); % Pitch rate                - rad/s  
% p           = data.Data(:,8); % Roll rate                 - rad/s
% r           = data.Data(:,9); % Yaw rate                  - rad/s
pitch       = deg2rad(pitch);% Pitch angle                  - rad
roll        = deg2rad(roll);% Roll angle                    - rad
heading_true= deg2rad(heading_true);% Heading               - rad
alpha       = deg2rad(alpha);% Angle of Attack              - rad
beta        = deg2rad(beta);% Side slip angle               - rad
latitude    = deg2rad(latitude);% Latitude                  - rad
longitude   = deg2rad(longitude);% Longitude                - rad  
altitude    = altitude * 0.3048;% Altitude                  - m 
% x           = data.Data(:,18);
% y           = data.Data(:,19);
% z           = data.Data(:,20);
% throttle_cmd = data.Data(:,21);% Throttle command         - %
% throttle_actual = data.Data(:,22);
eng_power   = eng_power * 745.7; % Engine power             - W
w_empty     = w_empty * 0.453592; 
w_payld     = w_payld * 0.453592; % Payload weight          - kg
w_fuel      = w_fuel * 0.453592; % Fuel weight              - kg 
% time        = data.Data(:,27); % Time from start          - seconds
% elevator    = data.Data(:,28); % Tail incidence           - % 
% aileron     = data.Data(:,29); %                          - %
% rudder      = data.Data(:,30); %                          - %


%% Testing Functions

[ ~, idx1 ] = min( abs( time-9359.83 ) );
[ ~, idx2 ] = min( abs( time-9534.34 ) );



V_test = vTrue_ktas(idx1:idx2);
t_test = time(idx1:idx2);

% Test Figure
figure;
plot(t_test, V_test, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks
[~,id] = findpeaks(V_test, 'MinPeakDistance', 20);
% displays peaks
plot(t_test(id), V_test(id), 'rx',MarkerSize=20)
xlabel('Time [seconds]')
ylabel('True Airspeed [ms^-^1]')
title('Time vs True Airspeed History TEST EXAMPLE')
grid on;
hold off;


[omega_damped, damping_ratio, omega_natural, lambda_1, lambda_2] = periodic_response(t_test, V_test);

lambda = non_periodic_response(t_test, V_test);
% 
% [mean_lambda, lambda] = non_periodic_response(time, vTrue_ktas)



%% SPPO (ossilitory)
% characterised by pitch

% test 1
[ ~, idx1 ] = min( abs( time-8784 ) );
[ ~, idx2 ] = min( abs( time-8807 ) );

Pitch_SPPO_1 = pitch(idx1:idx2);
t_SPPO_1 = time(idx1:idx2);

% figure;
% plot(t_SPPO_1, Pitch_SPPO_1, 'b', 'LineWidth', 1.5)
% hold on
% % finds Peaks
% [~,id] = findpeaks(Pitch_SPPO_1, 'MinPeakDistance', 5);
% % displays peaks
% plot(t_SPPO_1(id), Pitch_SPPO_1(id), 'rx',MarkerSize=20)
% xlabel('Time (seconds)')
% ylabel('Pitch [rad]')
% title('Time vs PitchHistory')
% grid on;
% hold off


SPPO_1 = periodic_response(t_SPPO_1, Pitch_SPPO_1);

% SPPO_1_lambda = non_periodic_response(t_SPPO_1, Pitch_SPPO_1);

%% test 2
[ ~, idx1 ] = min( abs( time-8812 ) );
[ ~, idx2 ] = min( abs( time-8814.9 ) );

Pitch_SPPO_2 = pitch(idx1:idx2);
t_SPPO_2 = time(idx1:idx2);

figure;
plot(t_SPPO_2, Pitch_SPPO_2, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks
[~,id1] = findpeaks(Pitch_SPPO_2, 'MinPeakDistance', 5);
[~,id2] = findpeaks(-Pitch_SPPO_2, 'MinPeakDistance', 5);


% displays peaks
plot(t_SPPO_2(id1), Pitch_SPPO_2(id1), 'rx',MarkerSize=20)
plot(t_SPPO_2(id2), Pitch_SPPO_2(id2), 'ro',MarkerSize=20)
xlabel('Time (seconds)')
ylabel('Pitch [rad]')
title('Time vs PitchHistory')
grid on;
hold off


SPPO_2 = periodic_response(t_SPPO_2, Pitch_SPPO_2);

% SPPO_2_lambda = non_periodic_response(t_SPPO_2, Pitch_SPPO_2);

%% test 3
[ ~, idx1 ] = min( abs( time-8835 ) );
[ ~, idx2 ] = min( abs( time-8843 ) );

Pitch_SPPO_3 = pitch(idx1:idx2);
t_SPPO_3 = time(idx1:idx2);

figure;
plot(t_SPPO_3, Pitch_SPPO_3, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks
[~,id] = findpeaks(Pitch_SPPO_3, 'MinPeakDistance', 5);
% displays peaks
plot(t_SPPO_3(id), Pitch_SPPO_3(id), 'rx',MarkerSize=20)
xlabel('Time (seconds)')
ylabel('Pitch [rad]')
title('Time vs PitchHistory')
grid on;
hold off


SPPO_3 = periodic_response(t_SPPO_3, Pitch_SPPO_3);

% SPPO_3_lambda = non_periodic_response(t_SPPO_3, Pitch_SPPO_3);

%% test 4
[ ~, idx1 ] = min( abs( time-8871 ) );
[ ~, idx2 ] = min( abs( time-8881 ) );

Pitch_SPPO_4 = pitch(idx1:idx2);
t_SPPO_4 = time(idx1:idx2);

figure;
plot(t_SPPO_4, Pitch_SPPO_4, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks
[~,id] = findpeaks(Pitch_SPPO_4, 'MinPeakDistance', 5);
% displays peaks
plot(t_SPPO_4(id), Pitch_SPPO_4(id), 'rx',MarkerSize=20)
xlabel('Time (seconds)')
ylabel('Pitch [rad]')
title('Time vs PitchHistory')
grid on;
hold off


SPPO_4 = periodic_response(t_SPPO_4, Pitch_SPPO_4);

% SPPO_4_lambda = non_periodic_response(t_SPPO_4, Pitch_SPPO_4);
% define data series


%% test 5
[ ~, idx1 ] = min( abs( time-11720 ) );
[ ~, idx2 ] = min( abs( time-11740 ) );

Pitch_SPPO_5 = pitch(idx1:idx2);
t_SPPO_5 = time(idx1:idx2);

% figure;
% plot(t_SPPO_5, Pitch_SPPO_5, 'b', 'LineWidth', 1.5)
% hold on
% % finds Peaks
% [~,id] = findpeaks(Pitch_SPPO_5, 'MinPeakDistance', 5);
% % displays peaks
% plot(t_SPPO_5(id), Pitch_SPPO_5(id), 'rx',MarkerSize=20)
% xlabel('Time (seconds)')
% ylabel('Pitch [rad]')
% title('Time vs PitchHistory')
% grid on;
% hold off


SPPO_5 = periodic_response(t_SPPO_5, Pitch_SPPO_5);

SPPO_5_lambda = non_periodic_response(t_SPPO_5, Pitch_SPPO_5);

%% test 6
[ ~, idx1 ] = min( abs( time-11750 ) );
[ ~, idx2 ] = min( abs( time-11760 ) );

Pitch_SPPO_6 = pitch(idx1:idx2);
t_SPPO_6 = time(idx1:idx2);

figure;
plot(t_SPPO_6, Pitch_SPPO_6, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks
[~,id] = findpeaks(Pitch_SPPO_6, 'MinPeakDistance', 5);
% displays peaks
plot(t_SPPO_6(id), Pitch_SPPO_6(id), 'rx',MarkerSize=20)
xlabel('Time (seconds)')
ylabel('Pitch [rad]')
title('Time vs PitchHistory')
grid on;
hold off


SPPO_6 = periodic_response(t_SPPO_6, Pitch_SPPO_6);

SPPO_6_lambda = non_periodic_response(t_SPPO_6, Pitch_SPPO_6);

%% test 7
[ ~, idx1 ] = min( abs( time-11780 ) );
[ ~, idx2 ] = min( abs( time-11790 ) );

Pitch_SPPO_7 = pitch(idx1:idx2);
t_SPPO_7 = time(idx1:idx2);

figure;
plot(t_SPPO_7, Pitch_SPPO_7, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks
[~,id] = findpeaks(Pitch_SPPO_7, 'MinPeakDistance', 5);
% displays peaks
plot(t_SPPO_7(id), Pitch_SPPO_7(id), 'rx',MarkerSize=20)
xlabel('Time (seconds)')
ylabel('Pitch [rad]')
title('Time vs PitchHistory')
grid on;
hold off


SPPO_7 = periodic_response(t_SPPO_7, Pitch_SPPO_7);

SPPO_7_lambda = non_periodic_response(t_SPPO_7, Pitch_SPPO_7);

%% test 8
[ ~, idx1 ] = min( abs( time-11800 ) );
[ ~, idx2 ] = min( abs( time-11820 ) );

Pitch_SPPO_8 = pitch(idx1:idx2);
t_SPPO_8 = time(idx1:idx2);

figure;
plot(t_SPPO_8, Pitch_SPPO_8, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks
[~,id] = findpeaks(Pitch_SPPO_8, 'MinPeakDistance', 5);
% displays peaks
plot(t_SPPO_8(id), Pitch_SPPO_8(id), 'rx',MarkerSize=20)
xlabel('Time (seconds)')
ylabel('Pitch [rad]')
title('Time vs PitchHistory')
grid on;
hold off


SPPO_8 = periodic_response(t_SPPO_8, Pitch_SPPO_8);

SPPO_8_lambda = non_periodic_response(t_SPPO_8, Pitch_SPPO_8);
%% Phugoid (ossilitory)

%% test 1
[ ~, idx1 ] = min( abs( time-9360 ) );
[ ~, idx2 ] = min( abs( time-9553 ) );

V_Phugoid_1 = vInd_kias(idx1:idx2);
t_Phugoid_1 = time(idx1:idx2);

figure;
plot(t_Phugoid_1, V_Phugoid_1, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks
[~,id1] = findpeaks(V_Phugoid_1, 'MinPeakDistance', 5);
[~,id2] = findpeaks(-V_Phugoid_1, 'MinPeakDistance', 5);


% displays peaks
plot(t_Phugoid_1(id1), V_Phugoid_1(id1), 'rx',MarkerSize=20)
plot(t_Phugoid_1(id2), V_Phugoid_1(id2), 'ro',MarkerSize=20)

xlabel('Time (seconds)')
ylabel('True Airspeed [ms^-1]')
title('Time vs True Airspeed')
grid on;
hold off


Phugoid_1 = periodic_response(t_Phugoid_1, V_Phugoid_1);

% Phugoid_1_lambda = non_periodic_response(t_Phugoid_1, V_Phugoid_1);


%% test 2
[ ~, idx1 ] = min( abs( time-9663 ) );
[ ~, idx2 ] = min( abs( time-9860 ) );

V_Phugoid_2 = vInd_kias(idx1:idx2);
t_Phugoid_2 = time(idx1:idx2);

figure;
plot(t_Phugoid_2, V_Phugoid_2, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks
[~,id1] = findpeaks(V_Phugoid_2, 'MinPeakDistance', 5);
[~,id2] = findpeaks(-V_Phugoid_2, 'MinPeakDistance', 5);


% displays peaks
plot(t_Phugoid_2(id1), V_Phugoid_2(id1), 'rx',MarkerSize=20)
plot(t_Phugoid_2(id2), V_Phugoid_2(id2), 'ro',MarkerSize=20)
xlabel('Time (seconds)')
ylabel('True Airspeed [ms^-1]')
title('Time vs True Airspeed')
grid on;
hold off


Phugoid_2 = periodic_response(t_Phugoid_2, V_Phugoid_2);

% Phugoid_2_lambda = non_periodic_response(t_Phugoid_2, V_Phugoid_2);

% define data series


%% Roll Test 

%% test 1
[ ~, idx1 ] = min( abs( time-9987 ) );
[ ~, idx2 ] = min( abs( time-10060 ) );

Roll_Roll_Subsidance_1 = roll(idx1:idx2);
t_Roll_Subsidance_1 = time(idx1:idx2);

figure;
plot(t_Roll_Subsidance_1, Roll_Roll_Subsidance_1, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks
[~,id1] = findpeaks(Roll_Roll_Subsidance_1, 'MinPeakDistance', 5);
[~,id2] = findpeaks(-Roll_Roll_Subsidance_1, 'MinPeakDistance', 5);
% displays peaks
plot(t_Roll_Subsidance_1(id2), Roll_Roll_Subsidance_1(id1), 'rx',MarkerSize=20)
plot(t_Roll_Subsidance_1(id2), Roll_Roll_Subsidance_1(id2), 'ro',MarkerSize=20)
xlabel('Time (seconds)')
ylabel('Roll [rad]')
title('Time vs Roll')
grid on;
hold off


Roll_Subsidance_1 = periodic_response(t_Roll_Subsidance_1, Roll_Roll_Subsidance_1);

Roll_Subsidance_1_lambda = non_periodic_response(t_Roll_Subsidance_1, Roll_Roll_Subsidance_1);

%% test 2
[ ~, idx1 ] = min( abs( time-10060 ) );
[ ~, idx2 ] = min( abs( time-10140 ) );

Roll_Roll_Subsidance_2 = roll(idx1:idx2);
t_Roll_Subsidance_2 = time(idx1:idx2);

figure;
plot(t_Roll_Subsidance_2, Roll_Roll_Subsidance_2, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks
[~,id1] = findpeaks(Roll_Roll_Subsidance_2, 'MinPeakDistance', 5);
[~,id2] = findpeaks(-Roll_Roll_Subsidance_2, 'MinPeakDistance', 5);
% displays peaks
plot(t_Roll_Subsidance_2(id1), Roll_Roll_Subsidance_2(id1), 'rx',MarkerSize=20)
plot(t_Roll_Subsidance_2(id2), Roll_Roll_Subsidance_2(id2), 'ro',MarkerSize=20)
xlabel('Time (seconds)')
ylabel('Roll [rad]')
title('Time vs Roll')
grid on;
hold off


Roll_Subsidance_2 = periodic_response(t_Roll_Subsidance_2, Roll_Roll_Subsidance_2);

Roll_Subsidance_2_lambda = non_periodic_response(t_Roll_Subsidance_2, Roll_Roll_Subsidance_2);


%% test 3
[ ~, idx1 ] = min( abs( time-10141 ) );
[ ~, idx2 ] = min( abs( time-10173 ) );

Roll_Roll_Subsidance_3 = roll(idx1:idx2);
t_Roll_Subsidance_3 = time(idx1:idx2);

figure;
plot(t_Roll_Subsidance_3, Roll_Roll_Subsidance_3, 'b', 'LineWidth', 1.5)
hold on
% finds Peaks
[~,id] = findpeaks(Roll_Roll_Subsidance_3, 'MinPeakDistance', 5);
% displays peaks
plot(t_Roll_Subsidance_3(id), Roll_Roll_Subsidance_3(id), 'rx',MarkerSize=20)
xlabel('Time (seconds)')
ylabel('Roll [rad]')
title('Time vs Roll')
grid on;
hold off


Roll_Subsidance_3 = periodic_response(t_Roll_Subsidance_3, Roll_Roll_Subsidance_3);

Roll_Subsidance_3_lambda = non_periodic_response(t_Roll_Subsidance_3, Roll_Roll_Subsidance_3);

% %% test 4
% [ ~, idx1 ] = min( abs( time-10290 ) );
% [ ~, idx2 ] = min( abs( time-10370 ) );
% 
% Roll_Roll_Subsidance_4 = roll(idx1:idx2);
% t_Roll_Subsidance_4 = time(idx1:idx2);
% 
% figure;
% plot(t_Roll_Subsidance_4, Roll_Roll_Subsidance_4, 'b', 'LineWidth', 1.5)
% hold on
% % finds Peaks
% [~,id] = findpeaks(Roll_Roll_Subsidance_4, 'MinPeakDistance', 5);
% % displays peaks
% plot(t_Roll_Subsidance_4(id), Roll_Roll_Subsidance_4(id), 'rx',MarkerSize=20)
% xlabel('Time (seconds)')
% ylabel('Roll [rad]')
% title('Time vs Roll')
% grid on;
% hold off
% 
% 
% Roll_Subsidance_4 = periodic_response(t_Roll_Subsidance_4, Roll_Roll_Subsidance_4);
% 
% Roll_Subsidance_4_lambda = non_periodic_response(t_Roll_Subsidance_4, Roll_Roll_Subsidance_4);
% 
% % define data series


%% Dutch Roll (ossilitory)

% define data series


%% Spiral Mode (non-ossilitory)

% define data series


%% Periodic Response Function

function [omega_damped, damping_ratio, omega_natural, lambda_1, lambda_2] = periodic_response(time, amplitude)

    % [pks,locs] = findpeaks(amplitude, 'MinPeakDistance', 20);


    [pks1, locs1] = findpeaks(amplitude, "MinPeakDistance", 5);
    [pks2, locs2] = findpeaks(-amplitude, "MinPeakDistance", 5);
    
    pks= sort([pks1;pks2], 'ascend');
    locs = sort([locs1;locs2], 'ascend');


    % findpeaks(amplitude);
    period = zeros( 1, length(locs) - 1 );
    delta = zeros( 1, length(locs) - 1 );

    for i = 1:( length(locs) - 1 )
        %calculates periods between peaks and throughs
        period(i) = 2*(time( locs(i+1) ) - time( locs(i) ));
        
        %clacluates logrithmic dectremnt
        delta(i) = log( 2*pks(i) / pks(i+1)); 
    end

    % delta = ln( amplitude(locs) / parameters(locs) );
    
    damping_ratio = mean(1 ./ sqrt( 1 + ( (2 * pi) ./ delta ) ));

    omega_damped = mean(2 * pi ./ period);

    omega_natural = mean(omega_damped ./ sqrt(1 - damping_ratio .^ 2));


    %calculates the eigenvalues of the state
    lambda_1 = ((-1) * omega_damped * damping_ratio)/(sqrt(1-damping_ratio^2)) - omega_damped*1i;
    lambda_2 = ((-1) * omega_damped * damping_ratio)/(sqrt(1-damping_ratio^2)) + omega_damped*1i;

end

%% Non-Periodic Response 

function [mean_lambda, lambda] = non_periodic_response(time, amplitude)

    [pks,locs] = findpeaks(amplitude);

    period = zeros( 1, length(locs) - 1 );
    
    lambda = zeros( 1, length(locs) - 1 );
    
    for i = 1:( length(locs) - 1 )
        %calculates periods between peaks
        period(i) = time( locs(i+1) ) - time( locs(i) );
        
        %calculates the eigenvalues of the state
        lambda(i) = (1 / period(i)) .* log(pks(i+1)./(pks(i)));

        % mean_lambda = mean(lambda);
    end
    mean_lambda = mean(lambda);

end

%% FFT for finding frequency of response


function [dominant_frequency] = calculateAverageFrequency(amplitude)

    % Load the data from the file
    % Extract the time column (first column)
    time = amplitude(:, 1);

    % Determine the number of vessels from the data
    num_vessels = (size(amplitude, 2) - 1) / 3; % Assuming y, v, h per vessel
    
    % Initialize a variable to store the sum of dominant frequencies
    dominant_frequency = zeros(num_vessels,1);
    
    % Loop through each vessel to calculate its dominant frequency
    for i = 1:num_vessels
    
        % Extract displacement data for the vessel
        displacement = amplitude(:, 3 * i - 1);
        
        % Perform Fourier Transform on the displacement data
        Y = fft(displacement);
        L = length(time);
        
        % Compute the two-sided spectrum and then the single-sided spectrum
        P2 = abs(Y / L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2 * P1(2:end-1);
        
        % Define the frequency domain
        f = (0:(L/2)) / max(time);
        
        % Find the frequency corresponding to the maximum amplitude
        [~, max_idx] = max(P1);
        dominant_frequency(i) = f(max_idx);
        
        % fprintf('The frequency of vessel %.1f is %.4f Hz\n',i ,dominant_frequency);
    end
    % Calculate the average frequency
    % Display the result
end
 


