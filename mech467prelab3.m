xPoints = [0 40 60];
yPoints = [0 30 30];
center = [90, 30];
R = 30;

Ti = 0.1*0.001;
A = 1000;
D= -A;
fs = 0;
fc = 200;
fe = 0;

L1 = sqrt(xPoints(2)^2+yPoints(2)^2);
L2 = sqrt((xPoints(3)-xPoints(2))^2+(yPoints(3)-yPoints(2))^2);
L3 = 2*pi*R;

T11 = (fc-fs)/A;
T31 = (fe-fc)/D;
T21 = L1/fc-fc/A;

N11 = T11/Ti;
N21 = T21/Ti;
N31 = T31/Ti;

k11 = 1:N11;
k21 = 1:N21;
k31 = 1:N31;

l11 = 0.5*A*(k11.*k11)*Ti^2;
l21 = fc*k21*Ti + l11(end);
l31 = 0.5*D*(k31.*k31)*Ti^2+fc*k31*Ti+l21(end);
f11 = A*k11*Ti;
f21 = fc*ones(1,fix(N21));
f31 = D*k31*Ti+fc;
a11 = A*ones(1,N11);
a21 = zeros(1,fix(N21));
a31 = D*ones(1,N31);


sinTheta = yPoints(2)/L1;
cosTheta = xPoints(2)/L1;
displacement = [l11 l21 l31];
velocity  = [f11 f21 f31];
acceleration = [a11 a21 a31];
time = [k11 k21+k11(end) k31+k21(end)+k11(end)];
xr=cosTheta*displacement;
yr = sinTheta*displacement;
xrprime = cosTheta*velocity;
yrprime = sinTheta*velocity;
xrprimeprime = cosTheta*acceleration;
yrprimeprime = sinTheta*acceleration;

% 





% 
% % first path
% 
% fcm = sqrt(2*A*D*L1/(D-A));
% T12 = fcm/A;
% T32 = fcm/A;
% 
% N12 = ceil(T12/Ti);
% N32 = ceil(T32/Ti);

% fcprime = 2*L1/(T12 +T32)
% Aprime = fcprime/T12;
% Dprime = -fcprime/T32;
% 
% k12 = 1:N12;
% k32 = 1:N32;
% 
% l12 = 0.5*Aprime*(k12.*k12)*Ti^2;
% l32 = 0.5*Dprime*(k32.*k32)*Ti^2+fcprime*k32*Ti+l12(end);
% f12  = Aprime*k12*Ti;
% f32  = Dprime*k32*Ti+fcprime;
% a12 = Aprime*ones(1,fix(N12));
% a32 = Dprime*ones(1,fix(N32));
% 
% displacement = [l12 l32];
% velocity = [f12 f32];
% acceleration = [a12 a32];
% time = [k12 k32+k12(end)];
% sinTheta = (yPoints(3)-yPoints(2))/L2;
% cosTheta = (xPoints(3)-xPoints(2))/L2;
% xr = [cosTheta*(l12) cosTheta*(l32)];
% yr = [sinTheta*(l12) sinTheta*(l32)];
% xrprime = [cosTheta*f12 cosTheta*f32];
% yrprime = [sinTheta*f12 sinTheta*f32];
% xrprimeprime = [cosTheta*a12 cosTheta*a32];
% yrprimeprime = [sinTheta*a12 sinTheta*a32];


%second path

fcm = sqrt(2*A*D*L2/(D-A));
T12 = fcm/A;
T32 = fcm/A;

N12 = ceil(T12/Ti);
N32 = T32/Ti;

fcprime = 2*L2/(T12 +T32)
Aprime = fcprime/T12;
Dprime = -fcprime/T32;

k12 = 1:N12;
k32 = 1:N32;

l12 = 0.5*Aprime*(k12.*k12)*Ti^2;
l32 = 0.5*Dprime*(k32.*k32)*Ti^2+fcprime*k32*Ti+l12(end);
f12  = Aprime*k12*Ti;
f32  = Dprime*k32*Ti+fcprime;
a12 = Aprime*ones(1,fix(N12));
a32 = Dprime*ones(1,fix(N32));

displacement = [displacement l12+displacement(end) l32+displacement(end)];
velocity = [velocity f12 f32];
acceleration = [acceleration a12 a32];
time = [time k12+time(end) k32+k12(end)+time(end)];
sinTheta = (yPoints(3)-yPoints(2))/L2;
cosTheta = (xPoints(3)-xPoints(2))/L2;
xr = [xr cosTheta*(l12)+xr(end) cosTheta*(l32)+xr(end)];
yr = [yr sinTheta*(l12)+yr(end) sinTheta*(l32)+yr(end)];
xrprime = [xrprime cosTheta*f12 cosTheta*f32];
yrprime = [yrprime sinTheta*f12 sinTheta*f32];
xrprimeprime = [xrprimeprime cosTheta*a12 cosTheta*a32];
yrprimeprime = [yrprimeprime sinTheta*a12 sinTheta*a32];


%third path

T13 = fc/A;
T33 = fc/A;
T23 = L3/fc-fc/A;

N13 = T13/Ti;
N23 = T23/Ti;
N33 = T33/Ti;

fcprime = 2*L3/(T13 + 2*T23 + T33);
Aprime = fcprime/T13;
Dprime = -Aprime;

xcyc = [90, 30];

k13 = 1:N13;
k23 = 1:N23;
k33 = 1:N33;

l13 = 0.5*Aprime*(k13.*k13)*Ti^2;
l23 = fcprime*k23*Ti+l13(end);
l33 = 0.5*Dprime*(k33.*k33)*Ti^2+fcprime*k33*Ti+l23(end);

f13 = Aprime*k13*Ti;
f23 = fcprime*ones(1,fix(N23));
f33 = Dprime*k33*Ti+fcprime;

a13 = Aprime*ones(1,N13);
a23 = zeros(1,fix(N23));
a33 = Dprime*ones(1,N33);


theta = pi + [l13 l23 l33]/R;
size([l13 l23 l33])
size(theta)
size([a13 a23 a33])

xr = [xr xcyc(1)+R*cos(theta)];
yr = [yr xcyc(2)+R*sin(theta)];
displacement = [displacement l13+displacement(end) l23+displacement(end) l33+displacement(end)];
time = [time k13+time(end) k23+k13(end)+time(end) k33+k13(end)+k23(end)+time(end)];
time = time*Ti;
velocity = [velocity f13 f23 f33];
acceleration = [acceleration a13 a23 a33];

xrprime = [xrprime [-f13 -f23 -f33].*sin(theta)];
yrprime = [yrprime [f13 f23 f33].*cos(theta)];
xrprimeprime = [xrprimeprime -[a13 a23 a33].*sin(theta)-1/R*[f13 f23 f33].^2.*cos(theta)];
yrprimeprime = [yrprimeprime [a13 a23 a33].*cos(theta)-1/R*[f13 f23 f33].^2.*sin(theta)];




% figure(4);
% 
% % First row of subplots
% subplot(2,3,1);
% plot(time, displacement);
% xlabel('Time (s)');
% ylabel('Displacement (mm)');
% title('Displacement vs Time');
% 
% subplot(2,3,2);
% plot(time, velocity);
% xlabel('Time (s)');
% ylabel('Velocity (mm/s)');
% title('Velocity vs Time');
% 
% subplot(2,3,3);
% plot(time, acceleration);
% xlabel('Time (s)');
% ylabel('Acceleration (mm/s^2)');
% title('Acceleration vs Time');
% 
% % Second row of subplots
% subplot(2,3,4);
% plot(time, xr);
% hold on;
% plot(time, yr);
% xlabel('Time (s)');
% ylabel('Position (mm)');
% title('Position vs Time');
% legend('x position', 'y position');
% 
% subplot(2,3,5);
% plot(time, xrprime);
% hold on;
% plot(time, yrprime);
% xlabel('Time (s)');
% ylabel('Velocity (mm/s)');
% title('Velocity vs Time');
% legend('x velocity', 'y velocity');
% 
% subplot(2,3,6);
% plot(time, xrprimeprime);
% hold on;
% plot(time, yrprimeprime);
% xlabel('Time (s)');
% ylabel('Acceleration (mm/s^2)');
% title('Acceleration vs Time');
% legend('x acceleration', 'y acceleration');
% 
% 
% 
% figure(5);
% plot(xr,yr);
% xlabel("xr (m)");
% ylabel("yr (m)");







% X-axis Lead Lag
wc = 2*pi*20; % Convert crossover frequency to rad/s
phase_margin = 60;    % Desired phase margin in degrees
% Plant transfer function (replace this with your actual plant transfer function)
numerator = [0.49*1.59];   % Numerator coefficients
denominator = [0.000436, 0.0094, 0];   % Denominator coefficients
plant_tfx = tf(numerator, denominator);



% Design lead-lag controller for 20Hz
controller_tf_x20 = leadlag(wc, phase_margin, 'X');

% Open-loop transfer function
open_loop_tf = controller_tf_x20 * plant_tfx;

% Plot Bode plot of the open-loop system
% figure(1)
% subplot(1,2,1);
% bode(open_loop_tf);
% title('Open Loop Bode Plot for X axis');
% hold on;
closed_loop_tf20 = feedback(open_loop_tf, 1);


wc = 2*pi*40; % Convert crossover frequency to rad/s
phase_margin = 60;    % Desired phase margin in degrees
controller_tf_x40 = leadlag(wc, phase_margin, 'X');

% Open-loop transfer function
open_loop_tf = controller_tf_x40 * plant_tfx;

closed_loop_tf40 = feedback(open_loop_tf, 1);





% subplot(1,2,2);
% bode(closed_loop_tf20);
% hold on;
% bode(closed_loop_tf40);
% bode(feedback(plant_tf,1));
% title('Closed Loop Step Bode Plot for X axis');
% legend('20Hz', '40Hz', "Uncompensated");

x20zeros = zero(closed_loop_tf20);
x20poles = pole(closed_loop_tf20);

x40zeros = zero(closed_loop_tf40);
x40poles = pole(closed_loop_tf40);


step_info = stepinfo(closed_loop_tf20);
overshoot = step_info.Overshoot;
rise_time = step_info.RiseTime;
bw = bandwidth(closed_loop_tf20);
x20step = [bw, overshoot, rise_time];

step_info = stepinfo(closed_loop_tf40);
overshoot = step_info.Overshoot;
rise_time = step_info.RiseTime;
bw = bandwidth(closed_loop_tf40);

x40step = [bw, overshoot, rise_time];


% figure(1);
% % lbwoutput = lsim(closed_loop_tf20, xr, time);
% % hbwoutput = lsim(closed_loop_tf40, xr, time);
% 
% plot(time, xr);
% hold on;
% plot(out.simout.time, out.simout.Data(:,2:5));
% 
% legend("X Target", "X Actual", "Y Target", "Y Actual");
% title("High BandWidth Target vs Simulated");
% xlabel("Time (s)");
% ylabel("Displacement (mm)");

% Y-axis Lead Lag
wc = 2*pi*20; % Convert crossover frequency to rad/s
phase_margin = 60;    % Desired phase margin in degrees

% Plant transfer function (replace this with your actual plant transfer function)
numerator = [0.49*1.59];   % Numerator coefficients
denominator = [0.0003, 0.0091, 0];   % Denominator coefficients
plant_tfy = tf(numerator, denominator);

% Design lead-lag controller for 20Hz
controller_tf_y20= leadlag(wc, phase_margin, 'Y');
% Open-loop transfer function
open_loop_tf = controller_tf_y20 * plant_tfy;
closed_loop_tf20 = feedback(open_loop_tf, 1);


wc = 2*pi*40; % Convert crossover frequency to rad/s
phase_margin = 60;    % Desired phase margin in degrees

% Design lead-lag controller for 40Hz
controller_tf_y40 = leadlag(wc, phase_margin, 'Y');
% Open-loop transfer function
open_loop_tf = controller_tf_y40 * plant_tfy;
closed_loop_tf40 = feedback(open_loop_tf, 1);


y20zeros = zero(closed_loop_tf20);
y20poles = pole(closed_loop_tf20);

y40zeros = zero(closed_loop_tf40);
y40poles = pole(closed_loop_tf40);

step_info = stepinfo(closed_loop_tf20);
overshoot = step_info.Overshoot;
rise_time = step_info.RiseTime;
bw = bandwidth(closed_loop_tf20);
y20step = [bw, overshoot, rise_time];
step_info = stepinfo(closed_loop_tf40);
overshoot = step_info.Overshoot;
rise_time = step_info.RiseTime;
bw = bandwidth(closed_loop_tf40);

y40step = [bw, overshoot, rise_time];

% figure(2);
% lbwoutput = lsim(closed_loop_tf20, yr, time);
% hbwoutput = lsim(closed_loop_tf40, yr, time);
% plot(time, yr);
% hold on;
% plot(time, lbwoutput);
% plot(time, hbwoutput);
% legend("Target", "Low Bandwidth Simulation", "High Bandwidth Simulation");
% title("Y axis Target vs Simulated");
% xlabel("Time (s)");
% ylabel("position (s)");




% x20poles
% x20zeros
% x40poles
% x40zeros
% 
% y20poles
% y20zeros
% y40poles
% y40zeros


%[[x20step];[x40step];y20step;y40step]

xcontroller  = controller_tf_x20;
ycontroller = controller_tf_y20;

% Open the Simulink model
open_system('prelab3')

% Set simulation parameters (if needed)
set_param('prelab3', 'StopTime', '2')

% Run the simulation
outExperiment = sim('prelab3');

% figure(1);
% plot(outExperiment.simout.Data(:,1), outExperiment.simout.Data(:,[6 7]));
% hold on;
% times = 0:Ti:47572*Ti;
% plot(times, Experimental(:,1)-Experimental(:,3));
% plot(times, Experimental(:,2)-Experimental(:,4));
% xlim([0 2]);
% 
% legend('X simulated error', 'Y simulated error', 'X actual error', 'Y actual error');
% xlabel("time (s)");
% ylabel("Error (mm)");








% figure(1);
% subplot(1,2,1);
% plot(outSlow.simout.Data(:,1), outSlow.simout.Data(:,[6]), 'DisplayName', 'Slow');
% hold on;
% plot(outFast.simout.Data(:,1), outFast.simout.Data(:,[6]),'DisplayName', 'Fast');
% 
% xlabel("Time (s)");
% ylabel("Error (mm)");
% legend();
% title("Tracking error at Low and High Feedrate in X");
% 
% subplot(1,2,2);
% plot(outSlow.simout.Data(:,1), outSlow.simout.Data(:,[7]),'DisplayName', 'Slow');
% hold on;
% plot(outFast.simout.Data(:,1), outFast.simout.Data(:,[7]),'DisplayName', 'Fast');
% legend();
% xlabel("Time (s)");
% ylabel("Error (mm)");
% 
% title("Tracking error at Low and High Feedrate in Y");


function Controller = leadlag(wc, phase_margin, type)
    % Design lead-lag controller
    if strcmpi(type, 'X')
        angle = phase_margin*pi/180;
        alpha = (1+sin(angle))/(1-sin(angle));
        Ts=1/wc/sqrt(alpha);
        Xnumerator = [0.49*1.59];   % Numerator coefficients
        Xdenominator = [0.000436, 0.0094, 0];   % Denominator coefficients
        Px = tf(Xnumerator, Xdenominator);
        Cx = tf([alpha*Ts 1],[Ts 1]);
        [LeadLagMagnitude, phase, frequency] = bode(Px*Cx, wc);
        Kp=1/LeadLagMagnitude;
        I =  tf([wc/10],[1 0]);
        Controller= Cx*Kp+I;
    elseif strcmpi(type, 'Y')
        angle = phase_margin*pi/180;
        alpha = (1+sin(angle))/(1-sin(angle));
        Ts=1/wc/sqrt(alpha);
        Ynumerator = [0.49*1.59];   % Numerator coefficients
        Ydenominator = [0.0003, 0.0091, 0];   % Denominator coefficients
        Py = tf(Ynumerator, Ydenominator);
        Cy = tf([alpha*Ts 1],[Ts 1]);
        [LeadLagMagnitude, phase, frequency] = bode(Py*Cy, wc);
        Kp=1/LeadLagMagnitude;
        I =  tf([wc/10],[1 0]);
        Controller = Cy*Kp+I;
    else
        error('Invalid type. Use ''Lead'' or ''Lag''.');
    end
    
end

