clear all; close all; clc;
% Here, we will see the complete process of retrieval of the original input signal 
% Through 16-QAM and will verify the same with theoretical values
% by plotting the Bit error rate(BER) vs Signal to Noise ratio(SNR) curve
% for naive coding.
M=16; %Modulation order
bb=log2(M); %Bits per symbol
loop=100; % No.of loops for frame data
N=100000; % Frame length (x_1 x_2 ... x_N)
x=randi([0,1],1,N); %Random input bits
nn=[];  % Initialization for storing normalized points
for i=1:4:length(x)
%*********************Naive encoding****************************************%    
    % Normalizing the points: sqrt(1/10) as
    % Average energy per constellation is 10
    if x(i)==0 && x(i+1)==0 && x(i+2)==0 && x(i+3)==0
        n=sqrt(1/10)*(-3+1j*3);
    elseif x(i)==0 && x(i+1)==1 && x(i+2)==0 && x(i+3)==0
        n=sqrt(1/10)*(-1+1j*3);  
    elseif x(i)==0 && x(i+1)==0 && x(i+2)==0 && x(i+3)==1
        n=sqrt(1/10)*(-3+1j*1);
    elseif x(i)==0 && x(i+1)==1 && x(i+2)==0 && x(i+3)==1
        n=sqrt(1/10)*(-1+1j*1);
        
    elseif x(i)==1 && x(i+1)==0 && x(i+2)==0 && x(i+3)==0
        n=sqrt(1/10)*(1+1j*3);   
    elseif x(i)==1 && x(i+1)==1 && x(i+2)==0 && x(i+3)==0
        n=sqrt(1/10)*(3+1j*3);
    elseif x(i)==1 && x(i+1)==0 && x(i+2)==0 && x(i+3)==1
        n=sqrt(1/10)*(1+1j*1);
    elseif x(i)==1 && x(i+1)==1 && x(i+2)==0 && x(i+3)==1
        n=sqrt(1/10)*(3+1j*1);
        
    elseif x(i)==0 && x(i+1)==0 && x(i+2)==1 && x(i+3)==0
        n=sqrt(1/10)*(-3+-1j*1);
    elseif x(i)==0 && x(i+1)==1 && x(i+2)==1 && x(i+3)==0
        n=sqrt(1/10)*(-1-1j*1);
    elseif x(i)==0 && x(i+1)==0 && x(i+2)==1 && x(i+3)==1
        n=sqrt(1/10)*(-3-1j*3);
    elseif x(i)==0 && x(i+1)==1 && x(i+2)==1 && x(i+3)==1
        n=sqrt(1/10)*(-1-1j*3);
        
    elseif x(i)==1 && x(i+1)==0 && x(i+2)==1 && x(i+3)==0
        n=sqrt(1/10)*(1-1j*1);
    elseif x(i)==1 && x(i+1)==1 && x(i+2)==1 && x(i+3)==0
        n=sqrt(1/10)*(3-1j*1);
    elseif x(i)==1 && x(i+1)==0 && x(i+2)==1 && x(i+3)==1
        n=sqrt(1/10)*(1-1j*3);
    elseif x(i)==1 && x(i+1)==1 && x(i+2)==1 && x(i+3)==1
        n=sqrt(1/10)*(3-1j*3);  
    end   
nn=[nn n]; % Storing the normalized points
end
BER_sim=[]; % Initialization for Storing BER of naive simulation
BER_th=[];  % Initialization for Storing BER of naive theoretical
for SNR_dB=0:0.5:16 % SNR in dB with step of 0.5
    SNR=10^(SNR_dB/10); % SNR in linear range
    noise=(1/sqrt(2))*[randn(1,length(nn))+1j*randn(1,length(nn))]; %AWGN noise
    sigma=sqrt(1/((log2(M))*SNR)); % deviation
    received=nn+sigma*noise;  % Channel
    ber=0; % Initiliasing ber for the given loop
    % Set of possible transmitted symbols
    trans=sqrt(1/10)*[(-3+1j*3) (-1+1j*3) (-3+1j*1) (-1+1j*1) (1+1j*3) (3+1j*3) (1+1j*1) (3+1j*1) (-3-1j*1) (-1-1j*1) (-3-1j*3) (-1-1j*3) (1-1j*1) (3-1j*1) (1-1j*3) (3-1j*3)];

    for lp=1:loop % Start looping of frame data
        %Minimum Hamming distance
        det=[];
        for m=1:length(received)       
            for l=1:length(trans)
                a=(real(received(m))-real(trans(l)))^2;
                b=(imag(received(m))-imag(trans(l)))^2;
                error(l)=sqrt(a+b);
            end
            iden=trans((error==min(error)));  % Finding the index for the minimum distance and identifying the transmitted symbol
            det=[det iden];
        end
        x_decoded=[];  % Initialisation for storing decoded bits
        for k=1:length(det)
%********************** Naive decoding ****************************************%
            if real(det(k))==sqrt(1/10)*-3 && imag(det(k))==sqrt(1/10)*3
                d=[0 0 0 0];
            elseif real(det(k))==sqrt(1/10)*-1 && imag(det(k))==sqrt(1/10)*3
                d=[0 1 0 0];
            elseif real(det(k))==sqrt(1/10)*-3 && imag(det(k))==sqrt(1/10)*1
                d=[0 0 0 1];
            elseif real(det(k))==sqrt(1/10)*-1 && imag(det(k))==sqrt(1/10)*1
                d=[0 1 0 1];
        
            elseif real(det(k))==sqrt(1/10)*1 && imag(det(k))==sqrt(1/10)*3
                d=[1 0 0 0];
            elseif real(det(k))==sqrt(1/10)*3 && imag(det(k))==sqrt(1/10)*3
                d=[1 1 0 0];
            elseif real(det(k))==sqrt(1/10)*1 && imag(det(k))==sqrt(1/10)*1
                d=[1 0 0 1];
            elseif real(det(k))==sqrt(1/10)*3 && imag(det(k))==sqrt(1/10)*1
                d=[1 1 0 1]; 
        
            elseif real(det(k))==sqrt(1/10)*-3 && imag(det(k))==sqrt(1/10)*-1
                d=[0 0 1 0];
            elseif real(det(k))==sqrt(1/10)*-1 && imag(det(k))==sqrt(1/10)*-1
                d=[0 1 1 0];
            elseif real(det(k))==sqrt(1/10)*-3 && imag(det(k))==sqrt(1/10)*-3
                d=[0 0 1 1];
            elseif real(det(k))==sqrt(1/10)*-1 && imag(det(k))==sqrt(1/10)*-3
                d=[0 1 1 1];
       
            elseif real(det(k))==sqrt(1/10)*1 && imag(det(k))==sqrt(1/10)*-1
                d=[1 0 1 0];
            elseif real(det(k))==sqrt(1/10)*3 && imag(det(k))==sqrt(1/10)*-1
                d=[1 1 1 0];
            elseif real(det(k))==sqrt(1/10)*1 && imag(det(k))==sqrt(1/10)*-3
                d=[1 0 1 1];
            elseif real(det(k))==sqrt(1/10)*3 && imag(det(k))==sqrt(1/10)*-3
                d=[1 1 1 1];    
            end   
        x_decoded=[x_decoded d]; % Storing decoded bits
        end
        rate=sum(x~= x_decoded)/N; % counter for checking mismatch of bits between input and output
        ber=ber+rate; % Adding the total BER till the end of the loop
    end   
    ber=ber/loop; % Average of BER from loop
    BER_sim=[BER_sim ber];                                                          % Storing BER of naive simulation for the range of SNR_dB
    BER_th=[BER_th (4/3)*(1/bb)*2*(1-sqrt(1/M))*erfc(sqrt((3*bb*SNR)/(2*(M-1))))]; % Storing BER of naive theoretical for the range of SNR_dB
end
SNR_dB=0:0.5:16;
semilogy(SNR_dB,BER_sim,'*m-',SNR_dB,BER_th,'g-'); % Plotting SNR vs BER for 16-QAM naive simulation & naive theoretical
axis([0 14 2*10^-5  1]);  % Setting the axis for the plots
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
grid on
legend('Simulated BER (Naive)','Theoretical BER (Naive)');
title('SNR vs BER plot for 16-QAM');

