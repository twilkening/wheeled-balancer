import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import statistics
import sys

def plot_serial_data(file_path):
    time = []
    gyro = []
    accelerometer = []

    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith("ms:"):
                continue  # skip lines that are not data
            try:
                parts = line.split(',')
                ms_str = parts[0].strip()
                time.append(float(ms_str.split(':')[1].strip()) / 1000.0)  # convert ms to seconds
                accelerometer.append(float(parts[1].split(':')[1].strip()))
                gyro.append(float(parts[2].split(':')[1].strip()))
            except ValueError:
                continue  # Skip lines that don't conform to expected format

    # compute average reading
    avg_gyro = sum(gyro) / len(gyro) if gyro else 0
    avg_accel = sum(accelerometer) / len(accelerometer) if accelerometer else 0
    print(f"Average Gyro: {avg_gyro:.2f}")
    print(f"Average Accelerometer: {avg_accel:.2f}")
    # compute standard deviation
    std_gyro = statistics.stdev(gyro) if len(gyro) > 1 else 0
    std_accel = statistics.stdev(accelerometer) if len(accelerometer) > 1 else 0
    print(f"Standard Deviation Gyro: {std_gyro:.3f}")
    print(f"Standard Deviation Accelerometer: {std_accel:.3f}")

    # Create separate figures for gyro and accelerometer
    # Gyro plot
    plt.figure(figsize=(12, 6))
    plt.subplot(2, 1, 1)
    plt.plot(time, gyro, marker='o', color='blue', linewidth=1, markersize=3)
    plt.xlabel('Time (s)')
    plt.ylabel('Gyro Value (deg/s)')
    plt.title('Gyro Readings Over Time')
    # plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    # plt.gca().xaxis.set_major_locator(mdates.SecondLocator(interval=1))  # Every 1 second
    # plt.xticks(rotation=45)
    plt.grid(True, alpha=0.3)
    
    # Accelerometer plot
    plt.subplot(2, 1, 2)
    plt.plot(time, accelerometer, marker='x', color='red', linewidth=1, markersize=3)
    plt.xlabel('Time (s)')
    plt.ylabel('Accelerometer Value (deg)')
    plt.title('Accelerometer Readings Over Time')
    # plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    # plt.gca().xaxis.set_major_locator(mdates.SecondLocator(interval=1))  # Every 1 second
    # plt.xticks(rotation=45)
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # Optional: Create separate windows
    # Uncomment the following if you prefer separate windows:
    
    # # Gyro in separate window
    # plt.figure(figsize=(10, 5))
    # plt.plot(ms, gyro, marker='o', color='blue', linewidth=1, markersize=3)
    # plt.xlabel('ms')
    # plt.ylabel('Gyro Value')
    # plt.title('Gyro Readings Over Time')
    # plt.xticks(rotation=45)
    # plt.grid(True, alpha=0.3)
    # plt.tight_layout()
    # plt.show()
    
    # # Accelerometer in separate window
    # plt.figure(figsize=(10, 5))
    # plt.plot(ms, accelerometer, marker='x', color='red', linewidth=1, markersize=3)
    # plt.xlabel('ms')
    # plt.ylabel('Accelerometer Value')
    # plt.title('Accelerometer Readings Over Time')
    # plt.xticks(rotation=45)
    # plt.grid(True, alpha=0.3)
    # plt.tight_layout()
    # plt.show()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        file_path = sys.argv[1]
    else:
        file_path = "log.csv"  # default value
    plot_serial_data(file_path)