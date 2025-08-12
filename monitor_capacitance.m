clc, clear, close all hidden;

% User configuration
deviceID = 'Dev1';        % Change if needed (use daq.getDevices)
aoChannel = 'ao0';
sampleRate = 1000;
rampTime = 1; % Ramp duration in seconds
initialTime = 1;
highHoldTime = 1;
lowHoldTime = 1;
maxVoltage = 4.0;
numActuations = 1;

% Use configuration for capacitance monitoring
aiVoltageChannel = 'ai4';
aiCurrentChannel = 'ai0';

% === Generate output signal ===
initialSamples = initialTime * sampleRate;
rampSamples = rampTime * sampleRate;
highHoldSamples = highHoldTime * sampleRate;
lowHoldSamples = lowHoldTime * sampleRate;
cycleSamples = 2*lowHoldSamples + 2*rampSamples + highHoldSamples;
totalSamples = initialSamples + numActuations*cycleSamples;
outputSignal = zeros(totalSamples,1);

% Initialize global variables
global inputData;
inputData = [];

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

% Create UI
fig = uifigure('Name', 'DAQ Voltage Ramp Control', 'Position', [100 100 700 500]);
ax1 = uiaxes(fig, 'Position', [50 280 600 180]);
title(ax1, 'Output Voltage');

ax2 = uiaxes(fig, 'Position', [50 60 600 180]);
title(ax2, 'Capacitance');

btnStart = uibutton(fig, 'push', 'Text', 'Start', ...
    'Position', [150 20 100 30], ...
    'ButtonPushedFcn', @(btn,event) startDAQ(deviceID, aoChannel, sampleRate, outputSignal, aiVoltageChannel, aiCurrentChannel));

btnSave = uibutton(fig, 'push', 'Text', 'Save CSV', ...
    'Position', [450 20 100 30], ...
    'ButtonPushedFcn', @(btn,event) saveCSV());


% Plot output signal preview
plot(ax1, t, outputSignal, 'b','LineWidth',1.5);
ax1.XLabel.String = 'Time (s)';
ax1.YLabel.String = 'Voltage (V)';
range = maxVoltage-min(outputSignal);
ylim(ax1, [min(outputSignal)-0.1*range, maxVoltage+0.1*range]);
xlim(ax1, [0,t(end)])


%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%
% Start DAQ function
function startDAQ(deviceID, aoChannel, sampleRate, outputSignal, aiVoltageChannel, aiCurrentChannel)
    dq = daq("ni");
    voltageOut = addoutput(dq, deviceID, aoChannel, 'Voltage');

    % Add input channel for recording sensor data (modify as needed)
    voltageIn = addinput(dq, "Dev1", aiVoltageChannel, "Voltage");
    voltageIn.TerminalConfig = "Differential";
    voltageIn.Name = "voltage_input";

    currentIn = addinput(dq, "Dev1", aiCurrentChannel, "Voltage");
    currentIn.TerminalConfig = "Differential";
    currentIn.Name = "current_input";

    dq.Rate = sampleRate;

    % Set the ScansAvailableFcn to store data
    dq.ScansAvailableFcn = @(src,evt) storeData(src, evt);

    % Set the 'ScansAvailableFcnCount' property
    dq.ScansAvailableFcnCount = floor(dq.Rate/10);

    % Zero out
    write(dq, [0]);

    outputSignal = double(outputSignal(:)); 
    outputSignal = [outputSignal];

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

    csvFileName = 'DAQ_data.csv';
    plotData(outputSignal, true, csvFileName);
end

% Function to store data
function storeData(src, event)
    global inputData;

    [eventData, eventTimestamps] = read(src, src.ScansAvailableFcnCount, "OutputFormat", "Matrix");

    if isempty(eventData)
        return; % If there's no data, simply return
    end

    newData = [eventTimestamps, eventData];
    inputData = [inputData; newData];

    % Report progress
    fprintf('Stored %d new data points. Total: %d\n', size(newData, 1), size(inputData, 1));
end

function plotData(outputSignal, saveData, csvFileName)
    global inputData;
    if isempty(inputData)
        uialert(fig, 'No data to plot yet.', 'Error');
        return;
    end
    
    % Extract voltage and current data
    % voltageData = inputData(:, 2) * 1e3; % kV -> V
    voltageData = outputSignal * 1e3; % kV -> V
    currentData = voltageToCurrent(inputData(:, 3)); % V -> A

    % Time vector
    timeVector = inputData(:, 1);

    % Accumulated charge
    accumulatedCharge = cumtrapz(timeVector, currentData);

    % Capacitance calculation
    capacitance = zeros(size(timeVector));
    valid_index = voltageData > 100; % Find valid indices
    capacitance(valid_index) = accumulatedCharge(valid_index) ./ voltageData(valid_index);

    % Plot capacitance
    figure;
    plot(timeVector, capacitance, 'g', 'LineWidth', 1.5);
    % % plot voltage
    % hold on;
    % plot(timeVector, voltageData, 'b', 'LineWidth', 1.5);
    % plot current
    % plot(timeVector, currentData, 'r', 'LineWidth', 1.5);
    % plot(timeVector, accumulatedCharge, 'r', 'LineWidth', 1.5);
    hold off;

    xlabel('Time (s)');
    ylabel('Capacitance (F)');
    title('Capacitance Over Time');
    grid on;

    if ~saveData
        % If not saving, just plot
        return;
    else
        saveCSV(voltageData, currentData, timeVector, capacitance, csvFileName);
    end
end

% Save CSV function
function saveCSV(voltageData, currentData, timeVector, capacitance, csvFileName)
    global inputData;
    if isempty(inputData)
        uialert(fig, 'No data to save yet.', 'Error');
        return;
    end
   
    %% Save data to CSV
    outputData = table(timeVector, voltageData, currentData, capacitance, ...
        'VariableNames', {'Time(s)', 'MeasuredVoltage(V)', 'Current(A)', 'Capacitance(F)'});

    writetable(outputData, csvFileName);
    disp(['Data saved to ' csvFileName]);
end

function convertedCurrentData = voltageToCurrent(currentData)
    % Convert voltage value of current monitor to current
    convertedCurrentData = currentData / 1 * 200 * 1e-6; % Convert from V to A
end

function capacitance = calculateCapacitance(voltageData, currentData, timeVector)
    % % Calculate accumulated charge
    % accumulatedCharge = cumtrapz(timeVector, currentData);

    % % Capacitance calculation
    % capacitance = zeros(size(timeVector));
    % valid_index = voltageData > 100; % Find valid indices
    % capacitance(valid_index) = accumulatedCharge(valid_index) ./ voltageData(valid_index);

    % Time step
    dt = diff(timeVector);

    % Voltage difference
    dV = diff(voltageData);

    % Charge difference (I * Δt)
    dQ = currentData(1:end-1) .* dt;

    % Capacitance calculation (ΔQ / ΔV)
    capacitance = zeros(size(voltageData));
    valid_idx = abs(dV) > 1e-3; % Exclude cases where voltage change is too small (threshold adjustable)
    % capacitance([false; valid_idx]) = dQ(valid_idx) ./ dV(valid_idx);
end