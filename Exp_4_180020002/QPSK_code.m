clear all; close all; clc;
% Here, we will see the complete process of retrieval of the original message signal 
% Through QPSK modulation and will verify the same with theoretical values
% by plotting both the Bit error rate(BER) vs Signal to Noise ratio(SNR) curves
% for both Grey coding and natural coding.
%--------------------------------------------------------------------------
loop=100; % No. of loops for frame data
M=10^4;  % Frame length (x_1 x_2 ... x_M)
SNR_dB=0:0.25:13; % SNR in dB with step of 0.25
SNR=10.^(SNR_dB./10); %SNR in linear range
BER_gray= zeros(size(SNR_dB)); % Initialising BER for gray mapping
BER_naive= zeros(size(SNR_dB)); % Initialising BER for natural mapping
Q_s=(1/sqrt(2))*[1+1j,-1+1j,-1-1j,1-1j]; %QPSK Symbols Constellation 
for snr= 1: length(SNR)   % start looping by SNR
    snr;
    for lp= 1: loop       % Start looping of frame data
        N_errors_gray=0;  % Error Counter for Gray 
        N_errors_Naive=0; % Error Counter for Naive
  
        for k=1:M
             
            x_input = unidrnd(length(Q_s));            % Random number generation
            T_s=Q_s(x_input);                          % Transmit symbols
            noise=(1/sqrt(2))*(randn(1)+1j*randn(1));  % Noise Generation
            R_s=sqrt(SNR(snr))*T_s+noise;              % receive symbols with noise
            
        
%****************Errors Computation for Gray Code (00 01 11 10) *************************%
            if sign(real(T_s))==sign(real(R_s)) && sign(imag(T_s))==sign(imag(R_s))
                N_errors_gray=N_errors_gray+0;
            elseif sign(real(T_s))==sign(real(R_s)) && sign(imag(T_s))~=sign(imag(R_s))
                N_errors_gray=N_errors_gray+1;
            elseif sign(real(T_s))~=sign(real(R_s)) && sign(imag(T_s))==sign(imag(R_s))
                N_errors_gray=N_errors_gray+1;
            else
                N_errors_gray=N_errors_gray+2;
            end
        
%****************Errors Computation for Naive Code (00 01 10 11)**************************%
            if sign(real(T_s))==sign(real(R_s)) && sign(imag(T_s))==sign(imag(R_s))
                N_errors_Naive=N_errors_Naive+0;
            elseif sign(real(T_s))~=sign(real(R_s)) && sign(imag(T_s))~=sign(imag(R_s))
                N_errors_Naive=N_errors_Naive+1;
            elseif sign(real(T_s))~=sign(real(R_s)) && sign(imag(T_s))==sign(imag(R_s))
                N_errors_Naive=N_errors_Naive+1;
            else
                N_errors_Naive=N_errors_Naive+2;
            end
        end    
        
% ********************* Bit Error Rate (BER) calulation ******************%
        BER_gray(snr)=BER_gray(snr)+N_errors_gray/(2*M);         % BER for Gray
        BER_naive(snr)=BER_naive(snr)+N_errors_Naive/(2*M);      % BER for Naive
    end        
    BER_gray(snr)= BER_gray(snr)/loop;                           % Average value over total no. of loops
    BER_naive(snr)= BER_naive(snr)/loop;                         % Average value over total no. of loops
                              
end 
% ******************* Plotting the simulation result ********************************************************% 
 BER_gray_theory=0.5*erfc(sqrt(SNR)/sqrt(2));         % BER for Gray thoery 
 BER_naive_theory=1.5*(0.5*erfc(sqrt(SNR)/sqrt(2)));  % BER for Naive thoery
 figure(1)
 semilogy(SNR_dB,BER_gray,'ob-','linewidth',2);        % Plot of BER for grey mapping simulation 
 hold on
 semilogy(SNR_dB,BER_gray_theory,'-r','linewidth',2); % Plot of BER for grey mapping theoretical
 hold on
 semilogy(SNR_dB,BER_naive,'*m-','linewidth',2);       % Plot of BER for natural mapping simulation
 hold on
 semilogy(SNR_dB,BER_naive_theory,'-k','linewidth',2);%Plot of BER for natural mapping theoretical
 axis([0 12 3*10^-5  1]);  
 xlabel( 'Signal to Noise Ratio (SNR)')
 ylabel( 'Bit Error Rate (BER)')
 title('Simulation QPSK transmission over noise');
 legend('Simulated BER (gray)', 'Theoretical BER (gray)','Simulated BER (Naive)', 'Theoretical BER (Naive)');
 grid on
% ************************* END ******************************************%




