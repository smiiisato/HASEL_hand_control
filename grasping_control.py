clc, clear, close all

% User configuration
deviceID = 'Dev1';        % Change if needed (use daq.getDevices)
aoChannel = 'ao0';
aiChannel = 'ai0';
sampleRate = 1000;
rampTime = 2;             % Ramp duration in seconds
highHoldTime = 4;
lowHoldTime = 2;
maxVoltage = 6.5;
numActuations = 5;

% === Generate output signal ===
initialSamples = highHoldTime * sampleRate;
rampSamples = rampTime * sampleRate;
highHoldSamples = highHoldTime * sampleRate;
lowHoldSamples = lowHoldTime * sampleRate;
cycleSamples = 2*lowHoldSamples + 2*rampSamples + highHoldSamples;
totalSamples = initialSamples + numActuations*cycleSamples;
outputSignal = nan(totalSamples,1);

% 1. Set inital pause
outputSignal(1:initialSamples) = zeros(initialSamples,1);

for i = 1:numActuations
    % 1. Initial low voltage
    secStart = initialSamples + (i-1)*cycleSamples;
    secEnd = secStart + lowHoldSamples - 1;
    outputSignal(secStart:secEnd) = zeros(lowHoldSamples,1);
    % 2. Ramp up
    secStart = secEnd + 1;
    secEnd = secStart + rampSamples - 1;
    rampUp = linspace(0, maxVoltage, rampSamples)';
    outputSignal(secStart:secEnd) = rampUp;
    % 3. Hold high voltage
    secStart = secEnd + 1;
    secEnd = secStart + highHoldSamples - 1;
    outputSignal(secStart:secEnd) = maxVoltage*ones(highHoldSamples,1);
    % 4. Ramp down
    secStart = secEnd + 1;
    secEnd = secStart + rampSamples - 1;
    rampDown = linspace(maxVoltage, 0, rampSamples)';
    outputSignal(secStart:secEnd) = rampDown;
    % 5. Low voltage pause
    secStart = secEnd + 1;
    secEnd = secStart + lowHoldSamples - 1;
    outputSignal(secStart:secEnd) = zeros(lowHoldSamples,1);
end

% time vector for plotting
t = (0:(totalSamples-1))' / sampleRate;
% Shared variables
sensorData = [];
timeStamps = [];
% Create UI
fig = uifigure('Name', 'DAQ Voltage Ramp Control', 'Position', [100 100 700 500]);
ax1 = uiaxes(fig, 'Position', [50 280 600 180]);
title(ax1, 'Output Voltage');

btnStart = uibutton(fig, 'push', 'Text', 'Start', ...
    'Position', [150 20 100 30], ...
    'ButtonPushedFcn', @(btn,event) startDAQ());

% Plot output signal preview
plot(ax1, t, outputSignal, 'b','LineWidth',1.5);
ax1.XLabel.String = 'Time (s)';
ax1.YLabel.String = 'Voltage (V)';
range = maxVoltage-min(outputSignal);
ylim(ax1, [min(outputSignal)-0.1*range, maxVoltage+0.1*range]);
xlim(ax1, [0,t(end)])

% Start DAQ function
function startDAQ()
    dq = daq("ni");
    addoutput(dq, deviceID, aoChannel, 'Voltage');
    dq.Rate = sampleRate;
    % Queue output signal
    preload(dq, outputSignal);
    start(dq, "Duration", seconds(length(outputSignal)/sampleRate));
    disp('DAQ started...');
end
