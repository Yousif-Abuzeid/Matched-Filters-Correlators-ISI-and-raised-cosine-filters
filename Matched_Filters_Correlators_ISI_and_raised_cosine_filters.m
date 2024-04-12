clear all
close all

%Generate Random data
num_bits_10=10;
samples_per_bit=5;
rand_data=generate_random_data(num_bits_10,samples_per_bit);
%Generate The Pulse Shaping Function
pulse = [5 4 3 2 1]/sqrt(55); 
%Convolve The Random Data With The Pulse
transmitter_out=conv(rand_data,pulse);
%Create a Matched Filter Using fliplr
matched_filter=fliplr(pulse);
%matched filter output
matched_filter_out=conv(transmitter_out,matched_filter);
%Create a Rect Filter
rect=[1 1 1 1 1]/sqrt(55);
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

%----requirement 1-(b)---- 
corr_out = correlator(transmitter_out,pulse,num_bits_10,samples_per_bit);
figure(2)
plot(matched_filter_out);
hold on;
plot(corr_out,'r');
title('Matched filter output Vs Correlator output'); 
xlabel('time');
ylabel('Amplitude');
hold off;




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

