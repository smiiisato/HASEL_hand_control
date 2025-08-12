% Close all figures, clear the command window, and clear variables
close all;
clc;
clear;

% Create a data acquisition object using a National Instruments device
dq = daq("ni");
flush(dq);

% Add a channel for outputting the commanded voltage
ch1 = addoutput(dq, "Dev1", "ao0", "Voltage");
ch1.TerminalConfig = "SingleEnded";

write(dq, [0.0]);

fprintf("it worked!\n");

close all;
% clc
clear 