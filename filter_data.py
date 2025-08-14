import pandas as pd
import numpy as np

MA = False
EMA = True

CSVFILE = 'ExperimentData/hand_grasp_6kV.csv'
CSVFILE2 = 'ExperimentData/hand_grasp_no_6kV.csv'
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
    # set capacitance as 0 if the voltage < 200
    data.loc[data['MeasuredVoltage(V)'] < 500, 'Capacitance(F)'] = 0

    # Apply moving average filter
    if MA:
        #filtered_capacitance = moving_average(data['Capacitance(F)'], window_size=5)
        capacitance = calculate_capacitance(data['MeasuredVoltage(V)'], data['Current(A)'], data['Time(s)'])
        filtered_capacitance = moving_average(capacitance, window_size=5)
    elif EMA:
        capacitance = calculate_capacitance(data['MeasuredVoltage(V)'], data['Current(A)'], data['Time(s)'])
        filtered_capacitance = exponential_moving_average(capacitance, alpha=0.01)
    return filtered_capacitance

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

    """ # --- Calculate accumulated charge (integral of I over time) ---
    accumulated_charge = cumtrapz(current_data, time_vector)

    # --- Capacitance calculation Q / V ---
    capacitance = np.zeros_like(time_vector, dtype=float)
    valid_index = voltage_data > 200
    capacitance[valid_index] = accumulated_charge[valid_index] / voltage_data[valid_index] """

    # --- Alternative method (commented out, ΔQ / ΔV) ---
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

def filter_data(file_path, file_path2=None):
    data = load_data(file_path)
    filtered_data = apply_filters(data)
    time = data['Time(s)']
    if file_path2:
        data2 = load_data(file_path2)
        filtered_data2 = apply_filters(data2)
        return time, filtered_data, filtered_data2
    return time, filtered_data, None

def plot_filtered_data(time, filtered_data, filtered_data2=None):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(12, 6))
    plt.plot(time, filtered_data, label=CSVFILE.split('/')[-1])
    if filtered_data2 is not None:
        plt.plot(time, filtered_data2, label=CSVFILE2.split('/')[-1])
    plt.title('Filtered Capacitance Data')
    plt.xlabel('Time')
    plt.ylabel('Capacitance')
    plt.legend()
    plt.grid()
    plt.show()

def plot_filtered_data_and_voltage(time, filtered_data, filtered_data2=None):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(12, 6))
    plt.plot(time, filtered_data, label=CSVFILE.split('/')[-1])
    if filtered_data2 is not None:
        plt.plot(time, filtered_data2, label=CSVFILE2.split('/')[-1])
    data = load_data(CSVFILE)
    plt.plot(time, data['MeasuredVoltage(V)']*1e-11, label=CSVFILE.split('/')[-1])
    plt.title('Filtered Capacitance Data')
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

    if CSVFILE2:
        data2 = load_data(CSVFILE2)
        plt.plot(time, data2[column_name], label=f'{column_name} - {CSVFILE2.split("/")[-1]}')

    plt.title(f'{column_name} Data')
    plt.xlabel('Time')
    plt.ylabel(column_name)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    time, filtered_data, filtered_data2 = filter_data(CSVFILE, CSVFILE2)
    #print(filtered_data)
    #plot_filtered_data(time, filtered_data, filtered_data2)
    plot_filtered_data_and_voltage(time, filtered_data, filtered_data2)
    #plot_data(time, "Current(A)")
    #plot_data(time, "MeasuredVoltage(V)")