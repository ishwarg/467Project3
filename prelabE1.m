subplot(2,2,1);
plot(xr, yr);
hold on;
xcontroller  = controller_tf_x20;
ycontroller = controller_tf_y20;

% Open the Simulink model
open_system('prelab3')

% Set simulation parameters (if needed)
set_param('prelab3', 'StopTime', '4')

% Run the simulation
out  = sim('prelab3');

plot(out.simout.Data(:,3), out.simout.Data(:,5));


ylabel("y (mm)");
xlabel("x (mm)");



xcontroller  = controller_tf_x40;
ycontroller = controller_tf_y40;

% Open the Simulink model
open_system('prelab3')

% Set simulation parameters (if needed)
set_param('prelab3', 'StopTime', '4')

% Run the simulation
out  = sim('prelab3');

plot(out.simout.Data(:,3), out.simout.Data(:,5));


xcontroller  = controller_tf_x40;
ycontroller = controller_tf_y20;

% Open the Simulink model
open_system('prelab3')

% Set simulation parameters (if needed)
set_param('prelab3', 'StopTime', '4')

% Run the simulation
out  = sim('prelab3');

plot(out.simout.Data(:,3), out.simout.Data(:,5));
legend("Target", "LBW", "HBW", "Mismatched");
xlim([18 22]);
ylim([12 16]);
title("R1 for different controllers");






subplot(2,2,2);
plot(xr, yr);
hold on;


xcontroller  = controller_tf_x20;
ycontroller = controller_tf_y20;

% Open the Simulink model
open_system('prelab3')

% Set simulation parameters (if needed)
set_param('prelab3', 'StopTime', '4')

% Run the simulation
out  = sim('prelab3');



plot(out.simout.Data(:,3), out.simout.Data(:,5));


ylabel("y (mm)");
xlabel("x (mm)");



xcontroller  = controller_tf_x40;
ycontroller = controller_tf_y40;

% Open the Simulink model
open_system('prelab3')

% Set simulation parameters (if needed)
set_param('prelab3', 'StopTime', '4')

% Run the simulation
out = sim('prelab3');

plot(out.simout.Data(:,3), out.simout.Data(:,5));


xcontroller  = controller_tf_x40;
ycontroller = controller_tf_y20;

% Open the Simulink model
open_system('prelab3')

% Set simulation parameters (if needed)
set_param('prelab3', 'StopTime', '4')

% Run the simulation
out = sim('prelab3');

plot(out.simout.Data(:,3), out.simout.Data(:,5));
legend("Target", "LBW", "HBW", "Mismatched");
xlim([38 42]);
ylim([28 32]);
title("R2 for different controllers");








subplot(2,2,3);
plot(xr, yr);
hold on;


xcontroller  = controller_tf_x20;
ycontroller = controller_tf_y20;

% Open the Simulink model
open_system('prelab3')

% Set simulation parameters (if needed)
set_param('prelab3', 'StopTime', '4')

% Run the simulation
out  = sim('prelab3');
plot(out.simout.Data(:,3), out.simout.Data(:,5));


ylabel("y (mm)");
xlabel("x (mm)");



xcontroller  = controller_tf_x40;
ycontroller = controller_tf_y40;

% Open the Simulink model
open_system('prelab3')

% Set simulation parameters (if needed)
set_param('prelab3', 'StopTime', '4')

% Run the simulation
out = sim('prelab3');

plot(out.simout.Data(:,3), out.simout.Data(:,5));


xcontroller  = controller_tf_x40;
ycontroller = controller_tf_y20;

% Open the Simulink model
open_system('prelab3')

% Set simulation parameters (if needed)
set_param('prelab3', 'StopTime', '4')

% Run the simulation
out = sim('prelab3');

plot(out.simout.Data(:,3), out.simout.Data(:,5));
legend("Target", "LBW", "HBW", "Mismatched");
xlim([58 62]);
ylim([28 32]);
title("R3 for different controllers");






subplot(2,2,4);
plot(xr, yr);
hold on;

xcontroller  = controller_tf_x20;
ycontroller = controller_tf_y20;

% Open the Simulink model
open_system('prelab3')

% Set simulation parameters (if needed)
set_param('prelab3', 'StopTime', '4')

% Run the simulation
out  = sim('prelab3');
plot(out.simout.Data(:,3), out.simout.Data(:,5));


ylabel("y (mm)");
xlabel("x (mm)");



xcontroller  = controller_tf_x40;
ycontroller = controller_tf_y40;

% Open the Simulink model
open_system('prelab3')

% Set simulation parameters (if needed)
set_param('prelab3', 'StopTime', '4')

% Run the simulation
out = sim('prelab3');

plot(out.simout.Data(:,3), out.simout.Data(:,5));


xcontroller  = controller_tf_x40;
ycontroller = controller_tf_y20;

% Open the Simulink model
open_system('prelab3')

% Set simulation parameters (if needed)
set_param('prelab3', 'StopTime', '4')

% Run the simulation
out = sim('prelab3');

plot(out.simout.Data(:,3), out.simout.Data(:,5));
legend("Target", "LBW", "HBW", "Mismatched");
xlim([110 114]);
ylim([48 52]);
title("R4 for different controllers")





figure(2)
plot(xr, yr);








