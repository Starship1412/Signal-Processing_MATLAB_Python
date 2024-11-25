from scipy.io import wavfile
import numpy as np
import sounddevice as sd
from playsound import playsound
# from IPython.display import Audio
import matplotlib.pyplot as plt

# Read in data
Fs, y = wavfile.read('pianonote.wav')
print("The sampling rate is %f" % Fs)  # print out the sampling rate of the audio file
x = y[:, 0]  # two channels are the same, so simply get the first channel
N = x.size  # number of samples
# By checking the data type of x, we see that x is in int16.
# That is each sample of x is represented by 16 bits (i.e. 2 bytes) in signed integer.
# Thus, the range of x is from -32768 to 32767.
# To further process x, we first need to convert it to the float32 format as follow
x = x/32767
plt.plot(x)
plt.figure()

# Compute the spectrum using FFT
X = np.fft.fft(x)
# Shift zero-frequency component to the centre of the spectrum
Xshifted = np.fft.fftshift(X)
Xshifted_mag = np.abs(Xshifted)  # magnitude spectrum
Xshifted_mag_dB = 20*np.log10(Xshifted_mag)  # magnitude spectrum in dB
Freq = np.fft.fftfreq(N, 1/Fs)  # x-axis frequencies
Freqshifted = np.fft.fftshift(Freq)

# Plot the magnitude spectrum in dB scale
plt.plot(Freqshifted, Xshifted_mag_dB)
plt.xlim(-1000, 1000)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude (dB)')
plt.show()