clear all; close all; clc;
% Here, we will see the complete process of retrieval of the original message signal 
% Through DQPSK modulation 
% by plotting the Bit error rate(BER) vs Signal to Noise ratio(SNR) curve
% for Non-Gray coding.
%--------------------------------------------------------------------------
loop=100; % No. of loops for frame data
N=4;
M=10000-1;  % Frame length (x_1 x_2 ... x_M), dont remove the -1 term
SNR_dB=0:1:12; % SNR in dB with step of 1
SNR=10.^(SNR_dB./10);
BER= zeros(1, length(SNR_dB)); % Initialising BER 
 
for dB= 1:length(SNR_dB)   % start looping by SNR
    for lp= 1: loop       % Start looping of frame data
%% *************Differential encoding*************************************%
        x_input= randi([0,1],1,M); % Random Input generation 
        dbit=1;
        x_diff=[]; % For storing encoded bits
        % Initial bit
        for i=1:length(x_input)
            x_diff=[x_diff dbit];
            dbit=not(xor(dbit,x_input(i))); % Calculating XNOR 
        end  
        x_diff=[x_diff dbit]; % Final encoded bits
%% ******************QPSK Modulation**************************************%
        %    Non-Gray Mapping
        %  Imaginary part
        %        ^
        %        |
        %  01*   |   *00
        %        |
        % <------+------> Real part
        %        |
        %  10*   |   *11
        %        |
        
        gg=[]; % % For storing mapped values from the bits 
        for i=1:2:length(x_diff)-1
            if x_diff(i)==0 && x_diff(i+1)==0        % 1st Quadrant
                g=cosd(45)+1j*sind(45);
            elseif x_diff(i)==0 && x_diff(i+1)==1    % 2nd Quadrant
                g=cosd(135)+1j*sind(135);
            elseif x_diff(i)==1 && x_diff(i+1)==0    % 3rd Quadrant
                g=cosd(225)+1j*sind(225);
            elseif x_diff(i)==1 && x_diff(i+1)==1    % 4th Quadrant
                g=cosd(315)+1j*sind(315);
            end
        gg=[gg g];
        end
%% ***************** Channel *************************************%    
        noise=(1/sqrt(2))*[randn(1,length(gg))+1j*randn(1,length(gg))]; %AWGN
        snr=10^(SNR_dB(dB)/10); % converting to linear scale
        sigma=sqrt(1/((log2(N)*snr))); %Standard deviation
        r=gg+sigma*noise;
%% ***************** QPSK Demodulation ***********************************%
        I=(imag(r)<0); % Inphase component
        Q=(not(sign(real(r))==sign(imag(r)))); % Perpendicular component
        x_out=[];
        for q=1:length(r)
            x_out=[x_out I(q) Q(q)]; % demodulated bits
        end
%% ***************** Differential decoding ********************************%
        for ii=1:length(x_out)-1
            if x_out(ii)==x_out(ii+1)  % 1 if Consecutive bits are equal
                det(ii)=1;
            else  
                det(ii)=0; % else 0
            end
        end
 %% ********************* Bit Error Rate (BER) calulation ******************%    
    
        ber=sum(x_input~=det)/M; % Computing ber
        BER(dB)= BER(dB) + ber; % add this ber to BER
    end
   
BER(dB)= BER(dB)/loop; % Average value over total no. of loops             
   
end 
%% ******************* Plotting the simulation result *********************%  
 f1 = figure(1);
 set(f1,'color',[1 1 1]);
 semilogy(SNR_dB,BER, 'b-o') % Plot for BER- Non-Gray simulation
 axis([0 10 10^-4  1.2]);  
 xlabel('Signal to Noise Ratio (SNR)')
 ylabel('Bit Error Rate (BER)')
 title('Simulation DQPSK transmission over noise');
 legend('BER simulation- Non-Gray')
 grid on;
%% ************************* END ******************************************%




