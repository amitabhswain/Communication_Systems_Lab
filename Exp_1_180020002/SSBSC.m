clear all; close all; clc;
% Here, we will see the complete process of retrieval of the original message signal
% Through SSB-SC modulation technique. 
%--------------------------------------------------------------------------
% Message signal
fs=1000; % Sampling frequency
t=0:1/fs:1; % Time
Am=2; fm=5; % Message Signal Amplitude, Frequency
message_signal= Am*cos(2*pi*fm*t); % Sinusoidal message signal
subplot(711); % Time domain Plot of Message Signal
plot(t,message_signal);
xlabel('Time');
ylabel('Amplitude');
title(['Sinusoidal message signal fm=', num2str(fm),'Hz']);
%--------------------------------------------------------------------------
% Carrier signals
% Carrier signal-1
Ac=3; fc=75; % Carrier Signal Amplitude, Frequency
carrier_signal_1= Ac*cos(2*pi*fc*t); %Sinusoidal carrier signal-1
subplot(712); %  Time domain Plot of Carrier Signal-1
plot(t,carrier_signal_1);
xlabel('Time');
ylabel('Amplitude');
title(['Sinusoidal Carrier Signal-1 fc=',num2str(fc),'Hz']);
% Carrier signal-2
Ac=3; fc=75; % Carrier Signal Amplitude, Frequency
carrier_signal_2= Ac*sin(2*pi*fc*t); %Sinusoidal carrier signal-2
subplot(713); %  Time domain Plot of Carrier Signal-2
plot(t,carrier_signal_2);
xlabel('Time');
ylabel('Amplitude');
title(['Sinusoidal Carrier Signal-2 fc=',num2str(fc),'Hz']);
%--------------------------------------------------------------------------
% Modulated signal
I=message_signal.*carrier_signal_1; % Inphase component
Q=imag(hilbert(message_signal)).*carrier_signal_2; % Orthogonal component
modulated_signal=I+Q; %Lower side band
%modulated_signal=I-Q; % Upper side band
subplot(714); %  Time domain Plot of Modulated Signal
plot(t,modulated_signal);
xlabel('Time');
ylabel('Amplitude');
title('SSBSC signal');
%--------------------------------------------------------------------------
% Frequency domain Plot of Modulated Signal
N=length(modulated_signal); % Number of DFT points
% Calculating the fft of modulated signal, then shifting the zero frequency component
% to center of the spectrum.
f_modulated_signal= fftshift(fft(modulated_signal,N));
f=fs*[-N/2:N/2-1]/N; % Frequency bins
subplot(715); % Plot
plot(f,abs(f_modulated_signal)); %Plotting the Magnitude values
xlabel('Absolute Frequency');
ylabel('DFT Values');
title('Frequency domain plot of SSBSC signal');
%--------------------------------------------------------------------------
% Demodulated Signal
% Multiply the modulated signal with 2*(carrier signal 1)
product_demodulator=2*modulated_signal.*carrier_signal_1; %Product demodulator output
% Low pass filtering
% Remove the 2*wc frequency components, then divide it by Ac^2 
% Butterworth filter design
% Mapping cutoff frequency range 0-fs/2 to 0-1[Normalisation]
% Setting cutoff frequency: 2*(fm/(fs/2))
% Setting filter order: 5
[num,den]=butter(5,4*fm/fs); % Calculating numerator, denominator of filter coefficients
demodulated_signal=filtfilt(num,den,product_demodulator); % Low pass filter with Zero phase shift
demodulated_signal=(demodulated_signal/(Ac*Ac));
subplot(716); %  Time domain Plot of demodulated signal 
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
subplot(717); % Plot
plot(f,abs(f_demodulated_signal));
xlabel('Absolute Frequency');
ylabel('DFT Values');
title('Frequency domain plot of demodulated signal');
%--------------------------------------------------------------------------

