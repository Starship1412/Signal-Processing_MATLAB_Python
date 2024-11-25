%% Task 1   Plot f(t)
clear all
clc
close all

t = -5:0.01:5;
a = 1;
ft = exp(-a*abs(t));

figure
plot(t,ft)
hold on
a = 2;
ft = exp(-a*abs(t));
plot(t,ft)
title('f(t)', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Time [s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)
legend('$a=1$', '$a=2$', 'Interpreter', 'Latex', 'FontSize', 10)
%% Task 2   Plot the CTFT of f(t)
clear all
clc
close all

Omega = -30:0.01:30;
a = 1;
Fw = 2*a./(a^2+Omega.^2); % expression of CTFT of the continuous time signal

figure
plot(Omega,Fw)
hold on
a = 2;
Fw = 2*a./(a^2+Omega.^2);
plot(Omega,Fw)
title('F($\Omega$)', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Frequency [rad/s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)
legend('$a=1$', '$a=2$', 'Interpreter', 'Latex', 'FontSize', 10)
%% Task 3   Sampling Period to Avoid Aliasing
clear all
clc
close all

t = -5:0.01:5;
a = 2;
ft = exp(-a*abs(t));
T = 1;
n = -5/T:5/T;
nT = n*T;
ft_ = exp(-a*abs(nT));

figure
plot(t,ft)
hold on
stem(nT,ft_)
title("f(t) " + "$\overline{f}$(t)", 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)
legend('f(t)', '$\overline{f}$(t)', 'Interpreter', 'Latex', 'FontSize', 10)

T = pi/30;
n = -5/T:5/T;
nT = n*T;
ft_ = exp(-a*abs(nT));

figure
plot(t,ft)
hold on
stem(nT,ft_)
title("f(t) " + "$\overline{f}$(t)", 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)
legend('f(t)', '$\overline{f}$(t)', 'Interpreter', 'Latex', 'FontSize', 10)
%% Task 4   first part   CTFT of the Sampled Signal
clear all
clc
close all

T = 1;
Omegas = 2*pi/T;
n = -5:5;
Omega = -2*Omegas:0.01:2*Omegas;
a = 2;
Fw_ = 0;

for k = 1:length(n)
    Fw_ = Fw_ + 1/T*2*a./(a^2+(Omega-n(k)*Omegas).^2); % expression of CTFT of the sampled signal
end

figure
plot(Omega,Fw_)
x_axis = [-2*Omegas -Omegas -1/2*Omegas 0 1/2*Omegas Omegas 2*Omegas];
xticks(x_axis)
title('$\overline{F}$($\Omega$)', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Frequency [rad/s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)
figure
plot(Omega,Fw_)
xlim([-1/2*Omegas 1/2*Omegas])
title('$\overline{F}$($\Omega$)', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Frequency [rad/s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)
%% Task 4   second part   CTFT of the Sampled Signal
clear all
clc
close all

T = pi/30;
Omegas = 2*pi/T;
n = -5:5;
Omega = -2*Omegas:0.01:2*Omegas;
a = 2;
Fw_ = 0;

for k = 1:length(n)
    Fw_ = Fw_ + 1/T*2*a./(a^2+(Omega-n(k)*Omegas).^2);
end

figure
plot(Omega,Fw_)
x_axis = [-2*Omegas -Omegas -1/2*Omegas 0 1/2*Omegas Omegas 2*Omegas];
xticks(x_axis)
title('$\overline{F}$($\Omega$)', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Frequency [rad/s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)
figure
plot(Omega,Fw_)
xlim([-1/2*Omegas 1/2*Omegas])
title('$\overline{F}$($\Omega$)', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Frequency [rad/s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)
%% Task 5   first part   DTFT of the Sampled Signal
clear all
clc
close all

a = 2;
T = 1;
Omega = -4*pi:0.01:4*pi;
num = 1 - exp(-2*a*T);
den = 1 - 2*cos(Omega)*exp(-a*T) + exp(-2*a*T);
F_w = num./den; % expression of DTFT of the sampled signal

figure
plot(Omega,F_w)
x_axis = [-4*pi -2*pi -pi 0 pi 2*pi 4*pi];
xticks(x_axis)
title('$\overline{F}$($\omega$)', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Frequency [rad/s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)
%% Task 5   second part   DTFT of the Sampled Signal
clear all
clc
close all

a = 2;
T = pi/30;
Omega = -4*pi:0.01:4*pi;
num = 1 - exp(-2*a*T);
den = 1 - 2*cos(Omega)*exp(-a*T) + exp(-2*a*T);
F_w = num./den;

figure
plot(Omega,F_w)
x_axis = [-4*pi -2*pi -pi 0 pi 2*pi 4*pi];
xticks(x_axis)
title('$\overline{F}$($\omega$)', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Frequency [rad/s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)
figure
plot(Omega,F_w)
xlim([-pi pi])
title('$\overline{F}$($\omega$)', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Frequency [rad/s]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)
%% Task 6   DTFT of the Truncated Sampled Signal
clear all
clc
close all

a = 2;
T = 0.1;
Omega = -pi:0.01:pi;

figure
hold on
N = 10;
n = -N:N;
ft_ = exp(-a*abs(n*T));
F1w = 0;
for n=-N:N
    F1w = F1w + ft_(n+N+1)*exp(-1j*n*Omega); % expression of DTFT of the truncated sampled signal
end
plot(Omega,abs(F1w))
N = 20;
n = -N:N;
ft_ = exp(-a*abs(n*T));
F2w = 0;
for n=-N:N
    F2w = F2w + ft_(n+N+1)*exp(-1j*n*Omega);
end
plot(Omega,abs(F2w))
num = 1 - exp(-2*a*T);
den = 1 - 2*cos(Omega)*exp(-a*T) + exp(-2*a*T);
F_w = num./den;
plot(Omega,F_w)
title('DTFT of the Truncated Sampled Signal compared with DTFT of the Sampled Signal', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)
legend('$\overline{F}_{1}$($\omega$)', '$\overline{F}_{2}$($\omega$)', '$\overline{F}$($\omega$)', 'Interpreter', 'Latex', 'FontSize', 10)
%% Task 7   DTFT and DFT
clear all
clc
close all

a = 2;
T = 0.1;
Omega = 0:0.01:2*pi;
N = 20;

figure
hold on
n = -N:N;
ft_ = exp(-a*abs(n*T));
F3 = zeros(2*N+1,1);
for m=0:2*N
    for n=-N:N
        F3(m+1) = F3(m+1) + ft_(n+N+1)*exp(-1j*2*pi*n.*m/(2*N)); % sum and assign values to each point of DFT
    end
end
w = 2*pi/(2*N)*((0:2*N)');
stem(w,abs(F3))
n = -N:N;
ft_ = exp(-a*abs(n*T));
F2w = 0;
for n=-N:N
    F2w = F2w + ft_(n+N+1)*exp(-1j*n*Omega); % expression of DTFT of the Truncated sampled signal
end
plot(Omega,abs(F2w))
title('DTFT and DFT', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('Time [ms]', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Amplitude', 'Interpreter', 'Latex', 'FontSize', 12)
legend('$F_{3}$', '$\overline{F}_{2}$($\omega$)', 'Interpreter', 'Latex', 'FontSize', 10)