%% Task 1   first part   Spectrum of a truncated cosine wave
clear all
clc
close all

A = 5;
f = 2e3;
fs = 8e3;
T = 1/fs;
N = 128;
td = (0:N-1)*T;
n=0:N-1;
x = A*cos(2*pi*f*td);
figure
stem(n,x)
title('x', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('n', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)

omega = 2*pi/N*(-N/2:N/2-1)*fs/(2*pi);
X = fft(x);
Xshifted = fftshift(X);
figure
stem(omega,abs(Xshifted)) % Draw the picture of absolute value of X considering that the zero-frequency component is shifted to the centre of the spectrum
title('X', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Frequency [rad/s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Magnitude', 'Interpreter', 'Latex', 'FontSize', 12)
%% Task 1   second part   Spectrum of a truncated cosine wave
clear all
clc
close all

A = 5;
f = 2e3;
fs = 8e3;
T = 1/fs;
N = 127;
td = (0:N-1)*T;
x = A*cos(2*pi*f*td);
figure
stem(td,x)
title('x', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Time [s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)

omega = 2*pi/N*(-N/2:N/2-1)*fs/(2*pi);
X = fft(x);
Xshifted = fftshift(X);
figure
stem(omega,abs(Xshifted)) % Draw the picture of absolute value of X considering that the zero-frequency component is shifted to the centre of the spectrum
title('X', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Frequency [rad/s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Magnitude', 'Interpreter', 'Latex', 'FontSize', 12)
%% Task 2   first part   Aliasing
clear all
clc
close all

A = 5;
f = 2e3;
fs = 3e3;
T = 1/fs;
N = 128;
td = (0:N-1)*T;
x = A*cos(2*pi*f*td);
figure
stem(td,x)
title('x', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Time [s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)

omega = 2*pi/N*(-N/2:N/2-1)*fs/(2*pi);
X = fft(x);
Xshifted = fftshift(X);
figure
stem(omega,abs(Xshifted)) % Draw the picture of absolute value of X considering that the zero-frequency component is shifted to the centre of the spectrum
title('X', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Frequency [rad/s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Magnitude', 'Interpreter', 'Latex', 'FontSize', 12)
%% Task 2   second part   Aliasing
clear all
clc
close all

A = 5;
f = 2e3;
fs = 3e3;
T = 1/fs;
N = 127;
td = (0:N-1)*T;
x = A*cos(2*pi*f*td);
figure
stem(td,x)
title('x', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Time [s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)

omega = 2*pi/N*(-N/2:N/2-1)*fs/(2*pi);
X = fft(x);
Xshifted = fftshift(X);
figure
stem(omega,abs(Xshifted)) % Draw the picture of absolute value of X considering that the zero-frequency component is shifted to the centre of the spectrum
title('X', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Frequency [rad/s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Magnitude', 'Interpreter', 'Latex', 'FontSize', 12)
%% Task 3   Zero Padding
clear all
clc
close all

A = 1;
f = 330.5;
fs = 1024;
T = 1/fs;
N = 2000; % 2023
td = (0:N-1)*T;
x = A*sin(2*pi*f*td);
omega = 2*pi/N*(-N/2:N/2-1);
X = fft(x);
Xshifted = fftshift(X);
figure
stem(2*pi/N*(-N/2:N/2-1), abs(Xshifted))
x_axis = [-pi -f/fs*2*pi 0 f/fs*2*pi pi];
xticks(x_axis)
title('$X$ \& $X_1$ \& $DTFT of x$', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Discrete-time frequency [rad/s] / Normalized Frequency [rad/s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Magnitude', 'Interpreter', 'Latex', 'FontSize', 12)

M = 2023;
x1 = [x.';zeros(M-N,1)];
omega_ = 2*pi/M*(-M/2:M/2-1);
X1 = fft(x1);
X1shifted = fftshift(X1);

[~,maxID] = max(abs(X1));
fmax_axis = (maxID-1)*2*pi/M % Compute the frequency associated with the peaks in the magnitude of X1

hold on
stem(omega_, abs(X1shifted))

x_DTFT = zeros(size(omega));
for k = 1:N
    x_DTFT = x_DTFT + x(k)*exp(-1i*omega*(k-1)); % Draw the DTFT of x according to the definition of DTFT
end
plot(omega, abs(x_DTFT), 'LineWidth', 1)
legend('X', 'X1', 'DTFT of x', 'Interpreter', 'Latex', 'FontSize', 10)