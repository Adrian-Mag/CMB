import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def STF(x):
    # Example STF function
    return np.exp(-((x - 25) / 9)**2)

# Generate some values for time
time = np.linspace(0, 50, 100)
dt = time[1] - time[0]

# Perform calculations and store the results
results = []
for t in time:
    result = [t, STF(t)]
    results.append(result)

# Compute the frequency spectrum
freq = np.fft.fftfreq(len(time), dt)
stf_freq = np.fft.fft(STF(time))
freq = np.fft.fftshift(freq)
stf_freq = np.fft.fftshift(stf_freq)

# Plot result
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=False)

# Plot STF function
ax1.plot(time, STF(time))
ax1.set_xlabel('Time [sec]')
ax1.set_ylabel('Amplitude')
ax1.set_xlim(time[0], time[-1])  # Set x-axis limits for ax1
ax1.grid(True)

# Plot frequency content
ax2.plot(freq, np.abs(stf_freq))
ax2.set_xlabel('Frequency [Hz]')
ax2.set_ylabel('Magnitude')
ax2.set_xlim(freq[0], freq[-1])  # Set x-axis limits for ax2
ax2.grid(True)
ax2.xaxis.set_major_locator(ticker.MaxNLocator(nbins=10))

plt.tight_layout()
plt.show()

# Ask if save
ans = input("Save the STF? (y/n): ")

if ans == "y":
    # Save results to a text file
    filename = "STF.txt"
    with open(filename, "w") as file:
        for result in results:
            file.write(f"{result[0]}\t{result[1]}\n")
