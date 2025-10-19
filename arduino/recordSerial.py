import serial
import time
from datetime import datetime

with serial.Serial('COM3', 57600, timeout=1) as ser:
    with open("log.csv", "w") as f:
        while True:
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]  # Include milliseconds
            line = ser.readline().decode(errors='ignore').strip()
            if line:
                # timestamped_line = f"{timestamp}, {line}"
                timestamped_line = f"{line}"
                print(timestamped_line)
                f.write(timestamped_line + "\n")
