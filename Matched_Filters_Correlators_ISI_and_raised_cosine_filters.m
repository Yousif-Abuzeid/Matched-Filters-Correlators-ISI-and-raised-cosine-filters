clear all
close all

%Generate Random data
num_bits_10=10;
samples_per_bit=5;
rand_data=generate_random_data(num_bits_10,samples_per_bit);
%Generate The Pulse Shaping Function
Bit_energy=55;
pulse = [5 4 3 2 1]/sqrt(Bit_energy); 
%Convolve The Random Data With The Pulse
transmitter_out=conv(rand_data,pulse);
%Create a Matched Filter Using fliplr
matched_filter=fliplr(pulse);
%matched filter output
matched_filter_out=conv(transmitter_out,matched_filter);
%Create a Rect Filter
rect_energy=5;
rect=[1 1 1 1 1]/sqrt(rect_energy);
%Rect Filter output
rect_filter_out=conv(transmitter_out,rect);

%plot both outcomes on the same figure
%----requirement 1-(a)----  
figure(1); 
subplot(2,1,1) 
plot(matched_filter_out);  
xlabel('time');
ylabel('Amplitude');
title('Matched filter output'); 
subplot(2,1,2)
xlabel('time');
ylabel('Amplitude');
plot(rect_filter_out, 'r');  
title('Rect filter output');

%plot both outcomes after sampling on the same figure
%matched filter
sampling_instants=5*[1:10];
sampled_matched_filter_out=matched_filter_out(5:5:50);
figure(2);
subplot(2,1,1)
stem(sampling_instants,sampled_matched_filter_out);
xlabel('time');
ylabel('Amplitude');
title('Matched filter output after sampling'); 
%rect filter
sampled_rect_filter_out=rect_filter_out(5:5:50);
subplot(2,1,2)
stem(sampling_instants,sampled_rect_filter_out,'r');
xlabel('time');
ylabel('Amplitude');
title('Rect filter output after sampling');

%----requirement 1-(b)---- 
corr_out = correlator(transmitter_out,pulse,num_bits_10,samples_per_bit);
figure(3)
plot(matched_filter_out);
hold on;
plot(corr_out,'r');
title('Matched filter output Vs Correlator output'); 
xlabel('time');
ylabel('Amplitude');
hold off;

%plot both outcomes after sampling 
figure(4);
%matched filter
stem(sampling_instants,sampled_matched_filter_out);
hold on;
%correlator
sampled_corr_out=corr_out(5:5:50);
stem(sampling_instants,sampled_corr_out,'r');
title('Matched filter output Vs Correlator output after sampling'); 
xlabel('time');
ylabel('Amplitude');
hold off;


%----2- Noise analysis: ---- 

%Generate Random data
num_bits_10K=10000;
samples_per_bit=5;
rand_data_10K=generate_random_data(num_bits_10K,samples_per_bit);
%Sample the signal every Ts (5 samples) to use it in BER caculations
sampled_transmitter_out_10K=rand_data_10K(1:5:50000);
%Convolve The Random Data With The Pulse
transmitter_out=conv(rand_data_10K,pulse);
%Create a Rect Filter
rect_energy=125;
rect_5=[5 5 5 5 5]/sqrt(rect_energy);
%Generate a unity variance, zero mean additive white Gaussian noise signal
%with the same size as transmitted signal.
noise=randn(size(transmitter_out));
%Loop on different values of SNR in dB
SNR_vector=[-2 -1 0 1 2 3 4 5];
normalized_energy_bit=1;
for i=1:length(SNR_vector)
    %resetting the noise
    noise=randn(size(transmitter_out));
    %calculate variance from SNR (SNR=Eb/No)
    No_vector(i)=normalized_energy_bit / 10^(SNR_vector(i)/10);
    variance_vector(i)=No_vector(i)/2; %variance=No/2
    %Scale the noise sequence to have variance = N0/2 by multiplying the sequence
    %by sqrt(N0/2). 
    noise=sqrt(variance_vector(i)) * noise; 
    %Add the noise to the transmitted sequence
    noisy_signal=transmitter_out + noise;
    %matched filter output
    matched_filter_out_10k=conv(noisy_signal,matched_filter);
    %Rect Filter output
    rect_filter_out_10K=conv(noisy_signal,rect_5);
    %Sample the matched filter output every Ts and estimate each bit  (5 samples)
    sampled_matched_filter_out_10k=estimate(samples_per_bit,num_bits_10K,matched_filter_out_10k);
    %Sample the rect filter output every Ts and estimate each bit  (5 samples)
    sampled_rect_filter_out_10k=estimate(samples_per_bit,num_bits_10K,rect_filter_out_10K);
    %Calculate the bit error rate for each SNR value for matched and rect
    %filters
    matched_error_counter=0;
    rect_error_counter=0;
    for c=1:length(sampled_matched_filter_out_10k)
        if sampled_transmitter_out_10K(c)~= sampled_matched_filter_out_10k(c)
             matched_error_counter=matched_error_counter + 1;
        end
        if sampled_transmitter_out_10K(c)~= sampled_rect_filter_out_10k(c)
             rect_error_counter=rect_error_counter + 1;
        end
    end
    BER_matched(i)=matched_error_counter/num_bits_10K;
    BER_rect(i)=rect_error_counter/num_bits_10K;
end
% calculate the theoritical BER
BER_theoritical=0.5 * erfc(sqrt(normalized_energy_bit ./ No_vector));
%plot matched filter BER vs theoritical
plot_BER(SNR_vector,BER_matched,BER_theoritical,'Matched filter output BER Vs theoritical BER');
%plot rect filter BER vs theoritical
plot_BER(SNR_vector,BER_rect,BER_theoritical,'rect filter output BER Vs theoritical BER');

% ------------ ISI and square root raised cosine -----------------%
% Parameters
num_bits = 100;
sample_per_bit = 5;
symbol_rate = 1; % sample per bit / 5 
data_train = generate_random_data(num_bits,sample_per_bit);
Rs = [0, 0, 1, 1];      % Rolloff factors for each case
delays = [2, 8, 2, 8];  % Filter delays for each case
num_cases = numel(Rs);

% Generate square root raised cosine filter coefficients for each case
for i = 1:num_cases
    R = Rs(i);
    delay = delays(i);
    filter =  rcosine(symbol_rate, sample_per_bit, 'sqrt', R, delay);
    tx_filtered_signal = conv(data_train, filter, 'same');
    rx_filtered_signal = conv(tx_filtered_signal, filter, 'same');
    % Plot eye diagram
    plot_eye_diagram(tx_filtered_signal, sample_per_bit, R, delay, 'Transmitter');
    plot_eye_diagram(rx_filtered_signal, sample_per_bit, R, delay, 'Receiver');
end



% Descripion :
% This Function Generates a Random sample Of Data
% Input : number of bits required 
%         number of samples per bit
% output : 
% Data: Random data in the form of +1 & -1 sampled per bit as given rate
%
function data = generate_random_data(num_bits,sample_per_bit)
%Generate Random data of ones & zeros
data = randi([0 1] , 1 , num_bits);
%Convert the Zeros to (-1)
data = (2*data) - 1; 
%upsample the Data using the sampling rate
data = upsample(data , sample_per_bit);
end

% Description :
% This Fucntion  recieves a signal and multiply it by a pulse and integrat 
%     along the resulting vector to calculate the correlation.
% Input: Recieved signal
%        our pulse   
%        number of bits required 
%        number of samples per bit
%output : the correlation between the pulse and recieved signal
function [corr_output] = correlator(Rx_signal , pulse , num_bits , sample_per_bit) 
%repeating the pulse signal to be properly multiplied. 
repeated_pulse=repmat(pulse,1,num_bits); 
%accumlator to apply the integration. 
accumlator =0; 
  
%for loop to integrate. 
for i = 1:num_bits*sample_per_bit 
     
    if mod((i-1),5)== 0 
        accumlator =0; 
    end 
   corr_output(i) = accumlator+ Rx_signal(i).*repeated_pulse(i); 
   accumlator= corr_output(i); 
  
end  
  
end

% Description :
% This Fucntion  Samples the matched filter output every Ts (5 samples)
%        and generates an array consisting of 10000 samples estimating the
%        value of each bit (1 or -1)
% Input: samples_per_bit
%        num_bits_10K   
%        matched_filter_out_10k 
%output: sampled_matched_filter_out_10k: the estimated array of bits
function [sampled_filter_out_10k]= estimate(samples_per_bit,num_bits_10K,filter_out_10k)
sampled_filter_out_10k=filter_out_10k(samples_per_bit:samples_per_bit:5*num_bits_10K);
for i=1:num_bits_10K
    if sampled_filter_out_10k(i)>=0
        sampled_filter_out_10k(i)=1;
    elseif sampled_filter_out_10k(i)<0
        sampled_filter_out_10k(i)=-1;
    end
end

end

% Description :
% This Fucntion plots the bit error rate for matched filter and rect filter
%        vs the theoritical BER
% Input: SNR_vector
%        BER_filter   
%        BER_theoritical 
%        str (recieves the title of the figure)
function plot_BER(SNR_vector,BER_filter,BER_theoritical,str)
figure;
semilogy(SNR_vector,BER_filter);
hold on;
semilogy(SNR_vector,BER_theoritical);
title(str); 
xlabel('Eb/No');
ylabel('BER');
hold off;
end
% Description:
% This function plots the eye diagram for a given filtered data with specified delay.
%
% Input:
% - filtered_data: The filtered data for which to plot the eye diagram.
% - delay: The delay used in the filtering process.
% - case_num: The number of the case (used for the title of the figure).
%
function plot_eye_diagram(filtered_data, sample_per_bit, R, delay, place)
    eyediagram(filtered_data, 2 * sample_per_bit);
    title(sprintf('%s Eye Diagram for Case: R=%d Delay=%d', place, R, delay));
    xlabel('Time');
    ylabel('Amplitude');
end