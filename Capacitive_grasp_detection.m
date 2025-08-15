clc, clear, close all hidden;

% User configuration
deviceID = 'Dev1';        % Change if needed (use daq.getDevices)
aoChannel = 'ao0';
sampleRate = 100000;
rampTime = 1; % Ramp duration in seconds
initialTime = 2;
highHoldTime = 0;
lowHoldTime = 0;
maxVoltage = 5.5;
numActuations = 1;

csvFileName = "hand_grasp_slow_.csv";
% csvFileName = "a.csv";

% Use configuration for capacitance monitoring
aiVoltageChannel = 'ai4';
aiCurrentChannel = 'ai0';

% === Generate output signal ===
global initialSamples rampSamples highHoldSamples lowHoldSamples cycleSamples totalSamples;
initialSamples = initialTime * sampleRate;
rampSamples = rampTime * sampleRate;
highHoldSamples = highHoldTime * sampleRate;
lowHoldSamples = lowHoldTime * sampleRate;
cycleSamples = 2*lowHoldSamples + 2*rampSamples + highHoldSamples;
totalSamples = initialSamples + numActuations*cycleSamples;

% Initialize global variables
global inputData;
inputData = [];

global outputSignal;
outputSignal = zeros(totalSamples,1);
outputSignal(1:initialSamples) = zeros(initialSamples,1);

global baselineCap;
baselineCap = 8.61 * 1e-8;

global capThreshold;
capThreshold = 3.5 * 1e-9;

global measuredAtMax;
measuredAtMax = false;

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

% setup DAQ
dq = setupDAQ(deviceID, aoChannel, sampleRate, outputSignal, aiVoltageChannel, aiCurrentChannel);

% time vector for plotting
t = (0:(totalSamples-1))' / sampleRate;

% Create UI
fig = uifigure('Name', 'DAQ Voltage Ramp Control', 'Position', [100 100 700 500]);
ax1 = uiaxes(fig, 'Position', [50 280 600 180]);
title(ax1, 'Output Voltage');

btnStart = uibutton(fig, 'push', 'Text', 'Start', ...
    'Position', [150 20 100 30], ...
    'ButtonPushedFcn', @(btn,event) startDAQ(dq, sampleRate, outputSignal, csvFileName));

% Plot output signal preview
plot(ax1, t, outputSignal, 'b','LineWidth',1.5);
ax1.XLabel.String = 'Time (s)';
ax1.YLabel.String = 'Voltage (V)';
range = maxVoltage-min(outputSignal);
ylim(ax1, [min(outputSignal)-0.1*range, maxVoltage+0.1*range]);
xlim(ax1, [0,t(end)])


%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%
% Start DAQ function
function startDAQ(dq, sampleRate, outputSignal, csvFileName)

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

    %plotData(outputSignal, true, csvFileName);
end

function dq = setupDAQ(deviceID, aoChannel, sampleRate, outputSignal, aiVoltageChannel, aiCurrentChannel)

    dq = daq("ni");
    voltageOut = addoutput(dq, "Dev1", aoChannel, 'Voltage');

    % Add input channel for recording sensor data (modify as needed)
    voltageIn = addinput(dq, "Dev1", aiVoltageChannel, "Voltage");
    voltageIn.TerminalConfig = "Differential";
    voltageIn.Name = "voltage_input";

    currentIn = addinput(dq, "Dev1", aiCurrentChannel, "Voltage");
    currentIn.TerminalConfig = "Differential";
    currentIn.Name = "current_input";

    dq.Rate = sampleRate;

    % Set the 'ScansAvailableFcnCount' property
    dq.ScansAvailableFcnCount = floor(dq.Rate/10);

    % Zero out
    write(dq, [0]);

    dq.ScansAvailableFcn = @(src, evt) storeData(src, evt);

    % Create a new figure for the capacitance real-time plot
    global hPlot capacitanceData timeData;
    capacitanceData = [];
    timeData = [];

    figure;
    hPlot = plot(NaN, NaN, 'LineWidth', 1.5);
    % add a horizontal line at baselineCap
    yline(baselineCap, 'r--', 'LineWidth', 1.5);
    xlabel('Sample');
    ylabel('Capacitance (F)');
    title('Grasp not detected', 'Color',[0, 0 ,1], 'FontSize', 20);
    grid on;

end

% Function to store data
function storeData(src, event)
    % global variables
    global inputData baselineCap capThreshold measuredAtMax outputSignal;

    % get data from DAQ
    [eventData, eventTimestamps] = read(src, src.ScansAvailableFcnCount, "OutputFormat", "Matrix");
    if isempty(eventData)
        return; % If there's no data, simply return
    end

    % Add new data
    newData = [eventTimestamps, eventData];
    inputData = [inputData; newData];
    fprintf('Stored %d new data points. Total: %d\n', size(newData, 1), size(inputData, 1));

    % Real-time plot of Capacitance
    real_time_plot_capacitance(inputData, newData, baselineCap, capThreshold, measuredAtMax);
end

function real_time_plot_capacitance(inputData, newData, baselineCap, capThreshold, measuredAtMax)
    % Real-time plot of Capacitance
    % This function plots the capacitance values in real-time
    % Add only new points to the plot when new data is available
    global hPlot capacitanceData timeData outputSignal;

    % === Set indices for voltage and current columns ===
    V = outputSignal(length(inputData)-length(newData):length(inputData), 1); % Voltage data
    I = newData(:,3);

    % === Capacitance estimation ===
    % Calculate capacitance using the formula C = I / (dV/dt) at each time this function is called
    dt = mean(diff(newData(:,1))); % Sampling interval [s]
    dVdt = [diff(V)/dt];
    %C = I ./ dVdt;
    C = emaFilterCapacitance(capacitanceData, I ./ dVdt, 0.1);
    %C(~isfinite(C)) = NaN; % Remove infinities/NaNs

    % Append data for plotting
    startIndex = length(capacitanceData) + 1;
    endIndex = startIndex + length(C) - 1;
    timeData = [timeData; newData(:,1)];
    capacitanceData = [capacitanceData; C];

    % Update plot
    set(hPlot, 'XData', timeData, 'YData', capacitanceData);
    if grasp_detection(inputData, V, C)
        title('Grasp detected', 'Color',[1, 0 ,0], 'FontSize', 20);
    else
        title('Grasp not detected', 'Color',[0, 0 ,1], 'FontSize', 20);
    end
end

function grasp_detected = grasp_detection(inputData, V, C)
    % Grasp detection logic
    global capThreshold measuredAtMax baselineCap initialSamples rampSamples;

    grasp_detected = false;
    if ~isempty(V) && length(inputData) >= initialSamples + rampSamples
        %Cmax = mean(C(V > max(V)*0.98));
        Cmax = max(C(V > max(V)*0.98));
        diffCap = baselineCap - Cmax;
        disp(['Cap change: ', num2str(diffCap*1e9), ' nF']);
        if diffCap > capThreshold
            grasp_detected = true;
        else
            disp('No grasp.');
            grasp_detected = false;
        end
        measuredAtMax = true;
    end
end

function filteredNewCap = emaFilterCapacitance(capacitanceData, newCap, alpha)
    % Exponential Moving Average filter
    filteredNewCap = zeros(size(newCap));
    if isempty(capacitanceData)
        filteredNewCap(1) = alpha * newCap(1); % Initialize with the first value
    else
        filteredNewCap(1) = (1 - alpha) * capacitanceData(end) + alpha * newCap(1);
    end

    for i = 2:length(capacitanceData)
        filteredNewCap(i) = (1 - alpha) * filteredNewCap(i-1) + alpha * newCap(i);
    end
end

function plotData(outputSignal, saveData, csvFileName)
    global inputData;
    if isempty(inputData)
        disp('No data to plot yet.');
        return;
    end
    
    % Extract voltage and current data
    voltageData = outputSignal * 1e3; % kV -> V
    currentData = voltageToCurrent(inputData(:, 3)); % V -> A

    % Time vector
    timeVector = inputData(:, 1);

    % Capacitance calculation
    capacitance = calculateCapacitance(voltageData, currentData, timeVector);

    % Plot capacitance
    figure;
    plot(timeVector, capacitance, 'g', 'LineWidth', 1.5);
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
        disp('No data to save yet.');
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
    % Time step
    dt = diff(timeVector);

    % Voltage difference
    dV = diff(voltageData);

    % Charge difference (I * Δt)
    dQ = currentData(1:end-1) .* dt;

    % Charge
    Q = cumtrapz(timeVector, currentData);

    % Capacitance calculation (ΔQ / ΔV)
    capacitance = zeros(size(voltageData));

    voltage_ramp_idx = abs(dV) > 1e-3; % Exclude cases where voltage change is too small (threshold adjustable)
    capacitance(voltage_ramp_idx) = dQ(voltage_ramp_idx) ./ dV(voltage_ramp_idx);

    voltage_step_idx = (abs(dV) <= 1e-3) & (voltageData(1:end-1) > 1000);
    capacitance(voltage_step_idx) = Q(voltage_step_idx) ./ voltageData(voltage_step_idx);
end