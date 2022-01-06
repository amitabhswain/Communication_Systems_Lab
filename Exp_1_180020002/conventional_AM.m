clear all; close all; clc;
% Here, we will see the complete process of retrieval of the original message signal
% Through Conventional AM. 
%--------------------------------------------------------------------------
% Message signal
fs=1000; % Sampling frequency
t=0:1/fs:1; % Time
Am=2; fm=5; % Message Signal Amplitude, Frequency
message_signal= Am*cos(2*pi*fm*t); % Sinusoidal message signal
subplot(611); % Time domain Plot of Message Signal
plot(t,message_signal);
xlabel('Time');
ylabel('Amplitude');
title(['Sinusoidal message signal fm=', num2str(fm),'Hz']);
%--------------------------------------------------------------------------
% Carrier signal
Ac=3; fc=75; % Carrier Signal Amplitude, Frequency
carrier_signal= Ac*cos(2*pi*fc*t); %Sinusoidal carrier signal 
subplot(612); %  Time domain Plot of Carrier Signal
plot(t,carrier_signal);
xlabel('Time');
ylabel('Amplitude');
title(['Sinusoidal Carrier Signal fc=',num2str(fc),'Hz']);
%--------------------------------------------------------------------------
% Modulated signal
% Representation of Conventional AM signal
modulated_signal=message_signal.*carrier_signal+carrier_signal;
subplot(613); %  Time domain Plot of Modulated Signal 
plot(t,modulated_signal);
xlabel('Time');
ylabel('Amplitude');
title('Conventional AM Signal');
%--------------------------------------------------------------------------
% Frequency domain Plot of Modulated Signal
N=length(modulated_signal); % Number of DFT points
% Calculating the fft of modulated signal, then shifting the zero frequency component
% to center of the spectrum.
f_modulated_signal= fftshift(fft(modulated_signal,N));
f=fs*[-N/2:N/2-1]/N; % Frequency bins
subplot(614); % Plot
plot(f,abs(f_modulated_signal)); %Plotting the Magnitude values
xlabel('Absolute Frequency');
ylabel('DFT Values');
title('Frequency domain plot of Conventional AM signal');
%--------------------------------------------------------------------------
% Demodulated Signal
% Multiply the modulated signal with 2*(carrier signal)
product_demodulator=2*modulated_signal.*carrier_signal; %Product demodulator output
% Low pass filtering
% Remove the 2*wc frequency components, then divide it by Ac^2 and subtract 1
% Butterworth filter design
% Mapping cutoff frequency range 0-fs/2 to 0-1[Normalisation]
% Setting cutoff frequency: 2*(fm/(fs/2))
% Setting filter order: 5
[num,den]=butter(5,4*fm/fs); % Calculating numerator, denominator of filter coefficients
demodulated_signal=filtfilt(num,den,product_demodulator); % Low pass filter with Zero phase shift
demodulated_signal=(demodulated_signal/(Ac*Ac))-1;
subplot(615); %  Time domain Plot of demodulated signal 
plot(t,demodulated_signal);
xlabel('Time');
ylabel('Amplitude');
title('Demodulated Signal');
%--------------------------------------------------------------------------
% Frequency domain Plot of demodulated signal
N=length(product_demodulator); %Number of DFT points
% Calculating the fft of demodulated signal, then shifting the zero frequency component
% to center of the spectrum.
f_demodulated_signal=fftshift(fft(demodulated_signal,N));
f=fs*[-N/2:N/2-1]/N; % Frequency bins
subplot(616); % Plot
plot(f,abs(f_demodulated_signal));
xlabel('Absolute Frequency');
ylabel('DFT Values');
title('Frequency domain plot of demodulated signal');
%--------------------------------------------------------------------------























