clc, clear, close all

% User configuration
deviceID = 'Dev1';        % Change if needed (use daq.getDevices)
aoChannel = 'ao0';
sampleRate = 1000;
rampTime = 1; % Ramp duration in seconds
initialTime = 1;
highHoldTime = 24;
lowHoldTime = 2;
maxVoltage = 5.5;
numActuations = 1;

% === Generate output signal ===
initialSamples = initialTime * sampleRate;
rampSamples = rampTime * sampleRate;
highHoldSamples = highHoldTime * sampleRate;
lowHoldSamples = lowHoldTime * sampleRate;
cycleSamples = 2*lowHoldSamples + 2*rampSamples + highHoldSamples;
totalSamples = initialSamples + numActuations*cycleSamples;
outputSignal = zeros(totalSamples,1);

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
    'ButtonPushedFcn', @(btn,event) startDAQ(deviceID, aoChannel, sampleRate, outputSignal));

% Plot output signal preview
plot(ax1, t, outputSignal, 'b','LineWidth',1.5);
ax1.XLabel.String = 'Time (s)';
ax1.YLabel.String = 'Voltage (V)';
range = maxVoltage-min(outputSignal);
ylim(ax1, [min(outputSignal)-0.1*range, maxVoltage+0.1*range]);
xlim(ax1, [0,t(end)])

% Start DAQ function
function startDAQ(deviceID, aoChannel, sampleRate, outputSignal)
    dq = daq("ni");
    addoutput(dq, deviceID, aoChannel, 'Voltage');
    dq.Rate = sampleRate;

    % Zero out
    write(dq, [0]);

    outputSignal = double(outputSignal(:)); 
    outputSignal = [outputSignal];
    disp(outputSignal);

    % Queue output signal
    preload(dq, outputSignal);
    start(dq, "Duration", seconds(length(outputSignal)/sampleRate));
    disp('DAQ started...');

    % Calculate expected duration and start acquisition
    expected_duration = length(outputSignal) / dq.Rate;
    disp(['Expected acquisition duration: ' num2str(expected_duration) ' seconds']);
    % Monitor acquisition progress
    start_time = tic;
    last_report_time = start_time;
    while toc(start_time) < expected_duration
        try
            pause(0.1);
            elapsed_time = toc(start_time);
            progress = (elapsed_time / expected_duration) * 100;
    
            % Report progress every second
            if toc(last_report_time) >= 1
                fprintf("Progress: %2.1f %%\n", progress);
                last_report_time = tic;
            end
        catch ME
            warning('Error in acquisition loop: %s', ME.message);
            break;
        end
    end
    
    stop(dq);
    pause(0.1);
    flush(dq);
    
    % Zero out
    write(dq, [0]);
end
