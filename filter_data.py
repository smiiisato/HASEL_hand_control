import pandas as pd
import numpy as np

MA = False
EMA = True

RAMP_TIME = 1  # Ramp duration in seconds
CSV_BASE = np.array([f'ExperimentData/baseline_ramp_{RAMP_TIME}s_{i}.csv' for i in range(1, 6)])
CSV_GRASP = np.array([f'ExperimentData/grasping_cube_ramp_{RAMP_TIME}s_{i}.csv' for i in range(1, 6)])
CSVFILE = np.concatenate((CSV_GRASP, CSV_BASE))

PLOT_CAPACITANCE_DIFF = True

CALCULATION_METHOD = 2
ALPHA = 0.0001
#CSVFILE2 = None

def moving_average(data, window_size):
    return data.rolling(window=window_size).mean()

def exponential_moving_average(values, alpha):
    ema = np.zeros_like(values, dtype=float)
    ema[0] = values[0] 
    for t in range(1, len(values)):
        ema[t] = alpha * values[t] + (1 - alpha) * ema[t-1]
    return ema

def apply_filters(data):
    # set capacitance as 0 if the voltage <= 0
    data.loc[data['MeasuredVoltage(V)'] <= 0, 'Capacitance(F)'] = 0

    # Apply moving average filter
    if MA:
        #filtered_capacitance = moving_average(data['Capacitance(F)'], window_size=5)
        capacitance = calculate_capacitance(data['MeasuredVoltage(V)'], data['Current(A)'], data['Time(s)'])
        filtered_capacitance = moving_average(capacitance, window_size=5)
    elif EMA:
        capacitance = calculate_capacitance(data['MeasuredVoltage(V)'], data['Current(A)'], data['Time(s)'])
        filtered_capacitance = exponential_moving_average(capacitance, alpha=ALPHA)
    return np.asarray(filtered_capacitance)

def calculate_capacitance(voltage_data, current_data, time_vector):
    """
    Calculate capacitance from voltage, current, and time data.

    Parameters
    ----------
    voltage_data : array_like
        Voltage values (V)
    current_data : array_like
        Current values (A)
    time_vector : array_like
        Time values (s)

    Returns
    -------
    capacitance : ndarray
        Calculated capacitance values
    """
    def cumtrapz(y, x):
        dx = np.diff(x)
        avg_y = (y[1:] + y[:-1]) / 2.0
        return np.concatenate(([0], np.cumsum(dx * avg_y)))

    voltage_data = np.asarray(voltage_data)
    current_data = np.asarray(current_data)
    time_vector = np.asarray(time_vector)

    # --- Calculate accumulated charge (integral of I over time) ---
    if CALCULATION_METHOD == 1:
        accumulated_charge = cumtrapz(current_data, time_vector)

        # --- Capacitance calculation Q / V ---
        capacitance = np.zeros_like(time_vector, dtype=float)
        valid_index = voltage_data > 200
        capacitance[valid_index] = accumulated_charge[valid_index] / voltage_data[valid_index]

    # --- Alternative method (commented out, ΔQ / ΔV) ---
    if CALCULATION_METHOD == 2:
        dt = np.diff(time_vector)
        dV = np.diff(voltage_data)
        dQ = current_data[:-1] * dt
        Q = cumtrapz(current_data, time_vector)

        capacitance = np.zeros_like(voltage_data, dtype=float)
        
        voltage_ramp_idx = np.abs(dV) > 1e-3
        capacitance[:-1][voltage_ramp_idx] = dQ[voltage_ramp_idx] / dV[voltage_ramp_idx]

        voltage_step_idx = (np.abs(dV) <= 1e-3) & (voltage_data[:-1] > 1000)
        capacitance[:-1][voltage_step_idx] = Q[:-1][voltage_step_idx] / voltage_data[:-1][voltage_step_idx]

    return capacitance


def load_data(file_path):
    return pd.read_csv(file_path)

def plot_filtered_data(time, CSVFILE):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(12, 6))
    for i, csvfile in enumerate(CSVFILE):
        filtered_data = apply_filters(load_data(csvfile))
        plt.plot(time, filtered_data, label=csvfile.split('/')[-1])
    plt.title('Filtered Capacitance Data')
    plt.xlabel('Time')
    plt.ylabel('Capacitance')
    plt.legend()
    plt.grid()
    plt.show()

def plot_filtered_data_and_voltage(CSVFILE):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(12, 6))
    for i, csvfile in enumerate(CSVFILE):
        filtered_data = apply_filters(load_data(csvfile))
        time = load_data(csvfile)['Time(s)']
        plt.plot(time, filtered_data, label=csvfile.split('/')[-1])
    data = load_data(CSVFILE[0])
    plt.plot(time, data['MeasuredVoltage(V)']*1e-11, label='Measured Voltage (V)', linestyle='--', color='gray')
    plt.title(f'calculation{CALCULATION_METHOD} - alpha={ALPHA}')
    plt.xlabel('Time')
    plt.ylabel('Capacitance')
    plt.legend()
    plt.grid()
    plt.show()

def plot_data(time, column_name):
    import matplotlib.pyplot as plt

    data = load_data(CSVFILE)
    plt.figure(figsize=(12, 6))
    plt.plot(time, data[column_name], label=CSVFILE.split('/')[-1])

    plt.title(f'{column_name} Data')
    plt.xlabel('Time')
    plt.ylabel(column_name)
    plt.legend()
    plt.show()

def calculate_capacitance_statistics(CSV_BASE, CSV_GRASP):
    """
    Calculate basic statistics for the capacitance data.

    Parameters
    ----------
    filtered_data : array_like
        Filtered capacitance values.

    Returns
    -------
    stats : dict
        Dictionary containing mean, median, and standard deviation.
    """
    data_base = [apply_filters(load_data(csvfile)) for csvfile in CSV_BASE]
    #print(f"Base data: {data_base}")
    data_grasp = [apply_filters(load_data(csvfile)) for csvfile in CSV_GRASP]
    #print(f"Grasp data: {data_grasp}")
    #print(f"Max Base data: {np.array([np.max(data) for data in data_base])}")
    base_max_value_mean = np.mean(np.array([np.max(data) for data in data_base]))
    base_max_value_median = np.median(np.array([np.max(data) for data in data_base]))
    base_max_value_max = np.max(np.array([np.max(data) for data in data_base]))
    base_max_value_min = np.min(np.array([np.max(data) for data in data_base]))
    grasp_max_value_mean = np.mean(np.array([np.max(data) for data in data_grasp]))
    grasp_max_value_median = np.median(np.array([np.max(data) for data in data_grasp]))
    grasp_max_value_max = np.max(np.array([np.max(data) for data in data_grasp]))
    grasp_max_value_min = np.min(np.array([np.max(data) for data in data_grasp]))

    stats = {
        "base_max_value_mean": base_max_value_mean,
        "base_max_value_median": base_max_value_median,
        "grasp_max_value_mean": grasp_max_value_mean,
        "grasp_max_value_median": grasp_max_value_median,
        "base_max_value_max": base_max_value_max,
        "base_max_value_min": base_max_value_min,
        "grasp_max_value_max": grasp_max_value_max,
        "grasp_max_value_min": grasp_max_value_min,
        "max_value_difference_mean": abs(grasp_max_value_mean - base_max_value_mean),
        "max_value_difference_median": abs(grasp_max_value_median - base_max_value_median),
        "max_value_difference_min": abs(grasp_max_value_max - base_max_value_min)
    }
    print(stats)

    if PLOT_CAPACITANCE_DIFF:
        import matplotlib.pyplot as plt
        time = load_data(CSVFILE[0])['Time(s)']
        plt.figure(figsize=(12, 6))
        base_filtered_data = np.array([apply_filters(load_data(csvfile)) for csvfile in CSV_BASE])
        grasp_filtered_data = np.array([apply_filters(load_data(csvfile)) for csvfile in CSV_GRASP])
        mean_base = np.mean(base_filtered_data, axis=0)
        mean_grasp = np.mean(grasp_filtered_data, axis=0)
        mean_diff = mean_base - mean_grasp
        plt.plot(time, mean_base, label='Base Mean')
        plt.plot(time, mean_grasp, label='Grasp Mean', linestyle='--')
        plt.plot(time, mean_diff, label='Mean Difference', linestyle=':')
        plt.title('Capacitance Data Comparison')
        plt.xlabel('Time (s)')
        plt.ylabel('Capacitance (F)')
        plt.legend()
        plt.grid()
        plt.show()
    return stats

if __name__ == "__main__":
    #print(filtered_data)
    plot_filtered_data_and_voltage(CSVFILE)
    calculate_capacitance_statistics(CSV_BASE, CSV_GRASP)