%% Assignment Requirements
% 1. The program should include the OFDM transmitter, 
% equivalent discrete-time channel, AWGN and OFDM demodulator. 
% 
% 2. For the equivalent discrete-time channel, 
% generate the channel filter coefficients h[n]μn=0 
% as independent and identically distributed zero-mean complex Gaussian random variables, 
% with variance 1/2 for real and imaginary parts.
% 
% 3. For the OFDM demodulator you can assume that the channel coefficients h[n]μn=0 are perfectly known.
% 
% 4. Test your MATLAB program with 64-QAM modulation on all data subcarriers and rc = 1/2 
% and obtain a plot of the bit error rate versus received Eb/N0.


%% Simulation
clear; clc; close; 

% arguments and constants
SNR = 1:1:40; % signal noise ratio (Eb/N0)
N_subcarrier = 64; 
scale = 6*2^16; % data scale

% convolutional encoder/decoder for rate 1/2 feedback
trellis = poly2trellis(5, [37 33], 37); 
tbdepth = 34; % Traceback depth for Viterbi decoder

% preallocation for speed
BER = zeros(1, length(SNR)); % bit error rate
convertor = [32*ones(1, 2*scale/6); 16*ones(1, 2*scale/6); 8*ones(1, 2*scale/6); ...
4*ones(1, 2*scale/6); 2*ones(1, 2*scale/6); ones(1, 2*scale/6)]; % bin2dec convertor

for snr = SNR
    %% Data generation
    % generate data bits
    bits = randi([0 1], scale, 1); 
    
    %% OFDM Transmitter
    % Encode
    enc_bits = convenc(bits, trellis); 
    
    % Reshape to [6, 2*scale/N_subcarrier] for the following bin2dec convertion
    enc_bits = reshape(enc_bits, [6, 2*scale/6]); 
    
    % Convert from binary to [0, 63] integers
    enc_bits = sum(enc_bits.*convertor); 
    
    % Modulate using 64-QAM
    complex_signals = qammod(enc_bits, 64); 
    
    % Reshape to N_subcarrier rows for the following 2^i-point ifft
    complex_signals = reshape(complex_signals, [N_subcarrier, 2*scale/(N_subcarrier*6)]); % reshape to a matrix who has 2^i row
    
    % Perform 2^i-point ifft operation
    complex_signals = ifft(complex_signals, N_subcarrier); 
    
    % Add cyclic prefix
    cp_complex_signals = zeros(80, 2*scale/(N_subcarrier*6)); 
    cp_complex_signals(17:end, :) = complex_signals; 
    cp_complex_signals(1:16, :) = complex_signals(49:64, :); 
    
    
    %% Equivalent discrete-time channel
    % Generate the channel filter coefficients h[n]μn=0
    % as independent and identically distributed zero-mean complex Gaussian random variables, 
    % with variance 1/2 for real and imaginary parts.
    h = 1/(sqrt(0.5*randn+0.5*randn*1i)); 
    channel_rayleigh = h*cp_complex_signals; 
    noise_gaussian = awgn(channel_rayleigh, snr-10*log10(6), 'measured'); % 6 bits every signal
    cp_complex_signals = h\noise_gaussian; 
    
    
    %% OFDM Receiver
    % Remove cyclic prefix
    complex_signals = cp_complex_signals(17:end, :); 
    
    % Perform 2^i-point fft operation
    complex_signals = fft(complex_signals, N_subcarrier); 
    
    % Reshape to one column
    complex_signals = reshape(complex_signals, [2*scale/6, 1]); 
    
    % Demodulate using 64-QAM
    enc_bits = qamdemod(complex_signals, 64); 
    
    % Convert from [0, 63] integers to binary
    enc_bits = dec2bin(enc_bits, 6)=='1'; 
    
    % Reshape to one column for the following decoding
    enc_bits = reshape(enc_bits', [scale*2, 1]); 
    
    % Decode
    dec_bits = vitdec(enc_bits, trellis, tbdepth, 'trunc', 'hard'); 
    
    
    %% BER Calculation
    [number, ratio] = biterr(dec_bits, bits); 
    BER(snr) = ratio; 
end

% plot
semilogy(SNR, BER); 
title("OFDM system bit error rate performances")
xlabel("SNR (Eb/N0)")
ylabel("BitErrorRate")





