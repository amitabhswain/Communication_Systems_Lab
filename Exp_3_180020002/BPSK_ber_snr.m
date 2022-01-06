clear all; close all; clc;
% Here, we will see the complete process of retrieval of the original message signal 
% Through BPSK modulation and will verify the same with theoretical values
% by plotting both the Bit error rate(BER) vs Signal to Noise ratio(SNR) curves
%--------------------------------------------------------------------------
loop=100; % No. of loops for frame data
M=100000;  % Frame length (x_1 x_2 ... x_M)
SNR_dB=0:0.1:12; % SNR in dB with step of 0.1
SNR=10.^(SNR_dB./10);
BER= zeros(1, length(SNR_dB)); % Initialising BER 
% ********************* Transmitter **************************************%
for dB= 1: length(SNR_dB) % start looping by SNR
    dB;
    for lp= 1: loop % Start looping of frame data
% ********************* BPSK signal generation ***************************%    
        x_input=sign(rand(M,1)- 0.5); % 1 for inphase and -1 for outphase
        N0=1./SNR(dB);
        sigma(dB)=sqrt(N0/2); % standard deviation
        noise=sigma(dB)*randn(M,1); % AWGN
% ********************* Channel ******************************************%
	    %y_channel=awgn(x_input,SNR_dB(dB)); % Additive White Gaussiann Noise (AWGN) 
        y_channel=x_input+noise;
    
% ********************* Receiver *****************************************% 
        y=y_channel;
        x_out= sign(real(y)); % Gives the output as 1, -1 depending on the sign of real part of y 
   
% ********************* Bit Error Rate (BER) calulation ******************%    
    
        [error, rate]= symerr(x_input, x_out); % Computes no. of symbol errors and symbol rate
        BER(dB)= BER(dB) + rate; % add this error rate to the BER
    end
   
BER(dB)= BER(dB)/loop; % Average value over total no. of loops
                              
   
end 
% ******************* Plotting the simulation result *********************%
   
 f1 = figure(1);
 set(f1,'color',[1 1 1]);
 semilogy(SNR_dB,BER, 'r') 
 hold on;
    
 BER_th=(1/2)*erfc(sqrt(SNR)); % theoretical calculation for BER
 semilogy(SNR_dB,BER_th,'g-o');   
 hold on;
 axis([0 12 0.000001  1.2]);  
 xlabel( 'Signal to Noise Ratio (SNR)')
 ylabel( 'Bit Error Rate (BER)')
 title('Simulation BPSK transmission over noise');
 legend('BER simulation','BER theoretical')
 grid on;
% ************************* END ******************************************%
    



