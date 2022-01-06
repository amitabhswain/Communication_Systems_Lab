clear all; close all; clc;
% Here, we will see the complete process of retrieval of the original message signal
% Through Frequency Modulation. 
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
% Frequency domain Plot of Message Signal
N=length(message_signal); % Number of DFT points
% Calculating the fft of message signal, then shifting the zero frequency component
% to center of the spectrum.
f_modulated_signal= fftshift(fft(message_signal,N));
f=fs*[-N/2:N/2-1]/N; % Frequency bins
subplot(612); % Plot
plot(f,abs(f_modulated_signal)); %Plotting the Magnitude values
xlabel('Absolute Frequency');
ylabel('DFT Values');
title('Frequency domain plot of message signal');
%--------------------------------------------------------------------------
% Carrier signal
Ac=3; fc=75; % Carrier Signal Amplitude, Frequency
%--------------------------------------------------------------------------
% Modulated signal
% Representation of FM signal
delta_f = 75;
kf = delta_f/fm; %frequency sensitivity of the modulator / modulation index
% Integration of m(tow) which is equal to sum_message
message_b = tril(ones(length(message_signal)));
message_c = message_signal.*message_b;
sum_message = sum(message_c,2);
dt = 1/fs;
% Expression for modulated signal
modulated_signal = Ac*cos(2*pi*fc*t+(2*pi*kf*sum_message').*dt);
subplot(613); %  Time domain Plot of Modulated Signal 
plot(t,modulated_signal);
xlabel('Time');
ylabel('Amplitude');
title('FM Signal');
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
title('Frequency domain plot of modulated signal');
%--------------------------------------------------------------------------
% Differentiating the modulated signal
differentiated_FM_signal = diff(modulated_signal); 
%Adding 0 to equalize the length of the vector
differentiated_FM_signal = [0,differentiated_FM_signal];
% Envelope detection
% Multipy the differentiated signal with 2 times the carrier signal
product_demodulator = 2*differentiated_FM_signal.*(Ac*(cos(2*pi*fc*t)));
% Low pass filtering
% Remove the 2*wc frequency components
% Butterworth filter design
% Mapping cutoff frequency range 0-fs/2 to 0-1[Normalisation]
% Setting cutoff frequency: 2*(fm/(fs/2))
% Setting filter order: 10
[num,den] = butter(10,4*fm/fs);
demodulated_signal = filtfilt(num,den,product_demodulator);
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
