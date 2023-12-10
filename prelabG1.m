plot(Experimental(:,3), Experimental(:,4));
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
legend("Experimental", "Simulated");
title("Experimental vs Simulated Toolpath");