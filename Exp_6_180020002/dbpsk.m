clear all; close all; clc;
% Here, we will see the complete process of retrieval of the original message signal 
% Through DBPSK modulation and will verify the same with theoretical values
% by plotting both the Bit error rate(BER) vs Signal to Noise ratio(SNR) curves
%--------------------------------------------------------------------------
loop=100; % No. of loops for frame data
M=10000;  % Frame length (x_1 x_2 ... x_M)
SNR_dB=0:1:12; % SNR in dB with step of 1
SNR=10.^(SNR_dB./10); % SNR in linear range
BER= zeros(1, length(SNR_dB)); % Initialising BER 
for dB= 1: length(SNR_dB) % start looping by SNR
    dB;
    for lp= 1: loop % Start looping of frame data
%% ********************* Differential encoding ***************************%
        x_input= randi([0,1],1,M); % Random Input generation
        dbit=1;
        x_diff=[]; % For storing encoded bits
        % Initial bit
        for i=1:length(x_input)
            x_diff=[x_diff dbit];
            dbit=not(xor(dbit,x_input(i))); % Calculating XNOR
        end  
        x_diff=[x_diff dbit]; % Final encoded bits
%% ********************* BPSK Modulation **********************************%        
        x_bpsk=sign(x_diff- 0.5); % 1 for inphase and -1 for outphase
%% ********************* Channel ******************************************%        
        N0=1./SNR(dB);
        sigma(dB)=sqrt(N0/2); % standard deviation
        noise=sigma(dB)*randn(length(x_bpsk),1)'; % AWGN
	    %y_channel=awgn(x_input,SNR_dB(dB)); % Additive White Gaussian Noise (AWGN) 
        y_channel=x_bpsk+noise;
%% ********************* Receiver *****************************************% 
        y=y_channel;
        x_rec= sign(real(y)); % Gives the output as 1, -1 depending on the sign of real part of y
%% ************** BPSK Demodulation **************************************%        
        for k=1:length(x_rec) % Mapping 1->1, -1->0
            if x_rec(k)==1
                x_rec(k)=1;
            else
                x_rec(k)=0;    
            end   
        end    
%% ********************* Differential decoding ***************************%
        for i=1:length(x_rec)-1
            if x_rec(i)==x_rec(i+1) % 1 if Consecutive bits are equal
                x_out(i)=1;
            else  
                x_out(i)=0; % else 0
            end
        end
   
%% ********************* Bit Error Rate(BER)calulation ******************%    
             
          ber=sum(x_input~=x_out)/M; % Computing ber
          BER(dB)= BER(dB) + ber; % add this ber to the BER
   end
   BER(dB)= BER(dB)/loop; % Average value over total no. of loops                  
   
end 
%% ******************* Plotting the simulation result *********************%
 f1 = figure(1);
 set(f1,'color',[1 1 1]);
 semilogy(SNR_dB,BER, 'r') % Plot for BER simulation
 hold on;
 BER_th=erfc(sqrt(SNR)); % theoretical calculation for BER 
 semilogy(SNR_dB,BER_th,'g-o'); % Plot for BER theoretical 
 hold on;
 axis([0 9 6*10^-5  1.2]); % Setting up the axis
 xlabel('Signal to Noise Ratio (SNR)')
 ylabel('Bit Error Rate (BER)')
 title('Simulation DBPSK transmission over noise');
 legend('BER simulation','BER theoretical')
 grid on;
%% ************************* END ******************************************%
    



