xPoints = [0 50 50];
yPoints = [0 50 50];
center = [90, 30];
R = 10;

Ti = 0.1*0.001;
A = 250;
D= -A;
fs = 0;
fc = 100;
fe = 0;

L1 = sqrt(xPoints(2)^2+yPoints(2)^2);
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


sinTheta = yPoints(2)/L1;
cosTheta = xPoints(2)/L1;
displacement = [l11 l21 l31];
time = [k11 k21+k11(end) k31+k21(end)+k11(end)];
xr=cosTheta*displacement;
yr = sinTheta*displacement;


%second path

T13 = fc/A;
T33 = fc/A;
T23 = L3/fc-fc/A;

N13 = T13/Ti;
N23 = T23/Ti;
N33 = T33/Ti;

fcprime = 2*L3/(T13 + 2*T23 + T33);
Aprime = fcprime/T13;
Dprime = -Aprime;

xcyc = [50+R 50];

k13 = 1:N13;
k23 = 1:N23;
k33 = 1:N33;

l13 = 0.5*Aprime*(k13.*k13)*Ti^2;
l23 = fcprime*k23*Ti+l13(end);
l33 = 0.5*Dprime*(k33.*k33)*Ti^2+fcprime*k33*Ti+l23(end);


theta = pi + [l13 l23 l33]/R;


xr = [xr xcyc(1)+R*cos(theta)];
yr = [yr xcyc(2)+R*sin(theta)];
displacement = [displacement l13+displacement(end) l23+displacement(end) l33+displacement(end)];
time = [time k13+time(end) k23+k13(end)+time(end) k33+k13(end)+k23(end)+time(end)];


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

xcyc = [50-R 50];

k13 = 1:N13;
k23 = 1:N23;
k33 = 1:N33;

l13 = 0.5*Aprime*(k13.*k13)*Ti^2;
l23 = fcprime*k23*Ti+l13(end);
l33 = 0.5*Dprime*(k33.*k33)*Ti^2+fcprime*k33*Ti+l23(end);


theta =[l13 l23 l33]/R;


xr = [xr xcyc(1)+R*cos(theta)];
yr = [yr xcyc(2)+R*sin(theta)];
displacement = [displacement l13+displacement(end) l23+displacement(end) l33+displacement(end)];
time = [time k13+time(end) k23+k13(end)+time(end) k33+k13(end)+k23(end)+time(end)];

plot(xr, yr);
xlabel("x (mm)")
ylabel("y (mm)")


time = time*Ti;

xrprime = diff(xr)/Ti;
yrprime = diff(yr)/Ti;
xrprimeprime = diff(xrprime)/Ti;
yrprimeprime = diff(yrprime)/Ti;

velocity = diff(displacement)/Ti;
acceleration = diff(velocity)/Ti;

traj.t = time;
traj.x = xr;
traj.y = yr;

save MyTraj traj

% Open the Simulink model
open_system('prelab3')

% Set simulation parameters (if needed)
set_param('prelab3', 'StopTime', '4')

% Run the simulation
out  = sim('prelab3');

figure(1);
plot(out.simout.Data(:,3), out.simout.Data(:,5));
hold on;
plot(Custom(:,1), Custom(:,2));
legend("Simulated", "Experimental");
xlabel("X (mm)");
ylabel("Y (mm)");
title("Experimental vs Simulated Toolpath of Custom Trajectory");


