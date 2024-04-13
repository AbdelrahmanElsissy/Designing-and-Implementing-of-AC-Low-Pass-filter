%% ------------------Exercise 4--------------------------
% Assignment 1
clc, close all; clear all;

u1 = 1; %V
R = 1e3; %Ohms
C = 100e-9; %F
f = logspace(1,5,100);
s = i*2*pi*f;
f_c = 1 / (R*C*2*pi);

%---Magnitude---
u2 = (u1) ./ ((2*pi*f.^2*R^2*C^2) + (3*R*C*s) + 1);
u3 = (u1) ./ (1 + i*2*pi*(R*2)*C*f);

figure('Position', [100 100 600 400]);
subplot(2,1,1)
loglog(f/1000, abs(u2), 'LineWidth', 2); hold on;
loglog(f/1000, abs(u3), 'LineWidth', 2);
xlabel("Frequency [kHz]"); ylabel("|H(f)|");
axis([20/1000 20 3e-2 2]);
legend("2nd Order LP Filter", "1st Order LP Filter")
H = line([0.0008 f_c/1000],[1 1]);
set(H, 'LineStyle', '-.', 'Color', 'r');
H = line([f_c/1000 150],[max(u2) 0.0008]);
set(H, 'LineStyle', '-.', 'Color', 'r');
text(f_c/1000,1.3,'f_c', 'FontSize', 8, 'Color', 'r');
text(0.1,0.5,'-20 dB/Dek', 'FontSize', 8, 'Color', 'r');

%----Phase----
subplot(2,1,2)
semilogx(f/1000, angle(u2)*180/pi, 'LineWidth', 2); hold on;
semilogx(f/1000, angle(u3)*180/pi, 'LineWidth', 2);
xlabel("Frequency (kHz)"); ylabel("Phase (H(f)) / Deg");
legend("2nd Order LP Filter", "1st Order LP Filter")
%--------------------------------------------------------------------------
%% ------------------------ assignment 2 -----------------------------------

clc, close all;
T = 3.5;  % Total signal duration in (s), more averages with longer signal
t_start = 0.4;                   % signal start: remove also transient part
t_stop=3.4;                      % signal end: remove the fading part

R = 1e3; %Ohms
C = 100e-9; %F
f = logspace(1,5,100);
s = i*2*pi*f;
f1_c = 1 / (R*C*2*pi);

%  Define parameters of your measurement
name='Test Transfer Function Measurement V0';  % name for output file
fs = 48000;             % sampling frequency, make sure value is correct!
N_freqs = 179;          % number of frequencies between f_start and f_stop 
f_start = 55;           % start frequency (Hz)
f_stop = 22000;         % stop frequency (Hz)
block_size = 2^12;      % block size (FFT window length)
ampl=0.01;              % select peak amplitude of output re full scale
%  Initialize audio interface and device reader
deviceReader = audioDeviceReader;
deviceReader.Driver = 'ASIO';

BufferSize = 1024;
saving = false;
aPR = audioPlayerRecorder('Device', 'ASIO4ALL v2', ...
    'SampleRate', fs,...
    'BitDepth', '16-bit integer',...
    'SupportVariableSize', true,...
    'BufferSize', BufferSize, ...
    'PlayerChannelMapping', [1,2], ...
    'RecorderChannelMapping', [1,2]);


% The number of output frequencies might be less than N_freqs, as the function removes redundant frequencies  
[t, sig, ms_indices] = generate_multisine(N_freqs, f_start, f_stop, fs, block_size, T);
sig=sig*ampl;  % scale signal
% ramp-up, ramp-down with hanning window
t_ramp = 200e-3;
n_ramp = floor(t_ramp*fs);
hann = hanning(2*n_ramp, 'periodic');
sig(1:n_ramp) = sig(1:n_ramp) .* hann(1:n_ramp);
sig(end-n_ramp+1:end) = sig(end-n_ramp+1:end) .* hann(n_ramp+1:end);


[playData, recData, N_underrun, N_overrun] = play_rec(aPR, sig);
release(aPR);

t_recData = (0:size(recData, 1)-1)'/fs;

% Process measurement data 
n_start = floor(t_start*fs);
n_stop = floor(t_stop*fs);
num_avg = floor((n_stop-n_start+1) /block_size); % number of blocks in data
n_stop = n_start + num_avg*block_size;          % align end to block border

rec = recData(n_start:n_stop-1,:);  % cut out data for analysis
t_rec = (0:length(rec)-1)/fs;

fprintf('Analyse data from T=%.2f to %.2f (averaging over %i blocks)\n',t_start,t_stop,num_avg);

% Average in time domain
rec_avg=mean(reshape(rec,block_size,num_avg,2),2);

% FFT of input and output
fft_rec = fft(rec_avg)/block_size;
fft_rec = 2*fft_rec(2:floor(block_size/2)+1,:); % single sided spectrum without DC
fft_freqs = fs*(1:block_size/2)'/block_size; % frequencies of spectrum
    
% Spectrum for frequencies, where energy was provided
ms_fft_freqs = fft_freqs(ms_indices,:); % select frequency vector
ms_fft_rec = fft_rec(ms_indices,:);     % select frequencies in spectrum

% Calculate transfer function
H=ms_fft_rec(:, 1) ./ ms_fft_rec(:, 2);  % Modified due to our connections
figure('Position', [100 100 800 400]);
subplot(2,1,1)
title('Transfer Function');
loglog(ms_fft_freqs, abs(H), 'b.'); 
L_fc = line([f1_c f1_c],[0.004 1]);
set(L_fc, 'LineStyle', '-.', 'Color', 'r');
text(f1_c+50 ,0.1,'f_c', 'FontSize', 8, 'Color', 'r');
ylabel('|H|');


subplot(2,1,2)
semilogx(ms_fft_freqs, unwrap(angle(H))*180/pi, 'b.');
xlabel('f [Hz]');
ylabel('Phase(H) [deg]');
L_fc_an= line([57 22000], [-75 -75]);
set(L_fc_an, 'LineStyle', '-.', 'Color', 'r');
text(f1_c-100 , -83, 'f_c', 'FontSize', 8, 'Color', 'r');
%--------------------------------------------------------------------------
%% ------------------------- Assignment 3-----------------------------------
clc, close all; clear all;

fs = 48000;
BufferSize = 1024;     % ASIO Buffer Size (Samples)
T = 5;
N = fs*T;
t = (0:N-1)'/fs;
sample_length = numel(t);
hann_win = fs * 0.1; 
window = ramp_filter(hann_win, sample_length);

noise = wgn(sample_length, 1, 1,1,1,'linear') / 12;
signal= noise .* window * sample_length/sum(window);

signal = signal(hann_win:sample_length-hann_win-1);
deviceReader = audioDeviceReader;
deviceReader.Driver = 'ASIO';
devices = getAudioDevices(deviceReader);


saving = false;

aPR = audioPlayerRecorder('Device', 'ASIO4ALL v2',...
    'SampleRate', fs,...
    'BitDepth', '16-bit integer',...
    'SupportVariableSize', true,...
    'BufferSize', BufferSize, ...
    'PlayerChannelMapping', [1,2], ...
    'RecorderChannelMapping', [1,2]);

[playData, recData, N_underrun,N_overrun] = play_rec(aPR, signal);


frequencies_ch2 = fs/numel(recData(:,2)) * (0:numel(recData(:,2))/2);
t_recorded = (0:(numel(recData(:,1))-1))'/fs;

figure('Position', [100 100 800 400]);
plot(t_recorded, recData(:,1), '-b', 'Linewidth', 1);  grid on; axis tight;
xlabel("Time [s]");
ylabel("Signal amplitude (V)");
figure('Position', [100 100 800 400]);
plot(t_recorded, recData(:,2), '-r', 'Linewidth', 1);  grid on; axis tight;
xlabel("Time [s]");
ylabel("Signal amplitude (V)");


% Mean and standard deviation of both channel recordings in time domain
mean_value1 = mean(recData(:,1))
mean_value2 = mean(recData(:,2))
std_dev1 = std(recData(:,1))
std_dev2 = std(recData(:,2))
%--------------------------------------------------------------------------
%% ----------------------- Assignment 4-------------------------------------
clc, close all; clear all;

fs = 48000;
BufferSize = 1024;     % ASIO Buffer Size (Samples)
T = 10; % Duration of the signal (s)
A_noise = 0.4; % Peak amplitude of white noise (V)
A_sine = 0.05; % Amplitude of sine signal (V)
f_sine = 375; % Frequency of sine signal (Hz)
N = fs*T;
t = (0:N-1)'/fs;
sample_length = numel(t);
hann_win = fs * 0.01; 
window = ramp_filter(hann_win, sample_length);
% Generate white noise with peak amplitude of 0.4 V
noise = A_noise * randn(size(t));

% Create sinusoidal signal with amplitude of 50 mV (0.05 V) and frequency of 375 Hz
sine_signal = A_sine * sin(2*pi*f_sine*t);

% Add sinusoidal signal to the noise
signal = noise + sine_signal;
f_signal= signal .* window * sample_length/sum(window);

% f_signal = f_signal(hann_win:sample_length-hann_win-1);
deviceReader = audioDeviceReader;
deviceReader.Driver = 'ASIO';
devices = getAudioDevices(deviceReader);


saving = false;

aPR = audioPlayerRecorder('Device', 'ASIO4ALL v2',...
    'SampleRate', fs,...
    'BitDepth', '16-bit integer',...
    'SupportVariableSize', true,...
    'BufferSize', BufferSize, ...
    'PlayerChannelMapping', [1,2], ...
    'RecorderChannelMapping', [1,2]);

[playData, recData, N_underrun,N_overrun] = play_rec(aPR, f_signal);


frequencies_ch2 = fs/numel(recData(:,2)) * (0:numel(recData(:,2))/2);
t_recorded = (0:(numel(recData(:,1))-1))'/fs;

figure('Position', [100 100 800 400]);
plot(t_recorded, recData(:,1), '-b', 'Linewidth', 1);  grid on; axis tight;
xlabel("Time [s]");
ylabel("Signal amplitude (V)");

figure('Position', [100 100 800 400]);
plot(t_recorded, recData(:,2), '-r', 'Linewidth', 1);  grid on; axis tight;
ylim([-1.2, 1.2]);
xlabel("Time [s]");
ylabel("Signal amplitude (V)");

% Define block size
block_size = 2^12; % 4096 samples

% Calculate the number of blocks
num_blocks = floor(numel(recData(:,2)) / block_size);

% Reshape the signal into a matrix with each column representing a block

output1 = recData(:,1);
output1_blocks = reshape(output1(1:num_blocks*block_size), block_size, num_blocks);
output2 = recData(:,2);
output2_blocks = reshape(output2(1:num_blocks*block_size), block_size, num_blocks);
% Calculate the average of each block
average_blocks1 = mean(output1_blocks');
average_blocks2 = mean(output2_blocks');

% Plot the average of blocks 1
figure('Position', [100 100 800 400]);
t_blocks = ((0:block_size-1) / fs) * num_blocks; 
% t_blocks = t_recorded / block_size;
plot(t_blocks, average_blocks1, 'color', 'blue'); hold on; axis tight;
xlabel('Time (s)');
ylabel('Average Amplitude (V)');
grid on;

% Plot the average of blocks 2
figure('Position', [100 100 800 400]);
plot(t_blocks, average_blocks2, 'color', 'red'); hold off; axis tight;
xlabel('Time (s)');
ylabel('Average Amplitude (V)');
grid on;

%mean and std
mean_value1 = mean(average_blocks1)
mean_value2 = mean(average_blocks2)
std_dev1 = std(average_blocks1) 
std_dev2 = std(average_blocks2)



% Compute the FFT for every block of the signal
Freq_bins = [1:1:118];
% signal_length = numel(recData(:,2));
for k = 1:num_blocks
    %average of complex spectra
    fft_blocks_double_sided1_raw = fft(output1_blocks(:,k));
    fft_blocks1_scale_raw = fft_blocks_double_sided1_raw/numel(output1_blocks(:,k));
    fft_blocks1_raw = fft_blocks1_scale_raw(1:numel(output1_blocks(:,k))/2+1);
    fft_blocks1_raw(2:end-1) = 2*fft_blocks1_raw(2:end-1);
    average_fft_blocks_raw(1,k) = mean(fft_blocks1_raw);
    %average for absolute values
    fft_blocks_double_sided1 = fft(output1_blocks(:,k));
    fft_blocks1_scale = abs(fft_blocks_double_sided1/numel(output1_blocks(:,k)));
    fft_blocks1 = fft_blocks1_scale(1:numel(output1_blocks(:,k))/2+1);
    fft_blocks1(2:end-1) = 2*fft_blocks1(2:end-1);
    average_fft_blocks(1,k) = mean(fft_blocks1);


end

figure('Position', [100 100 800 400]);
plot(Freq_bins, abs(average_fft_blocks_raw));
xlabel('Number of Bins');
ylabel('Magnitude (V)');
grid on; axis tight;

Freq_bins = Freq_bins(5:(numel(Freq_bins)-2));
average_fft_blocks =  average_fft_blocks(5:(numel(average_fft_blocks)-2));

figure('Position', [100 100 800 400]);
plot(Freq_bins, average_fft_blocks);
xlabel('Number of Bins');
ylabel('Magnitude (V)');
grid on; axis tight;

% Define block size
block_size1 = 2^9; % 4096 samples
% Calculate the number of blocks
num_blocks1 = floor(numel(recData(:,2)) / block_size1);
% Reshape the signal into a matrix with each column representing a block
output11_blocks = reshape(output1(1:num_blocks1*block_size1), block_size1, num_blocks1);
Freq_bins1 = [1:1:num_blocks1];
for k = 1:num_blocks1
    fft_blocks_double_sided11 = fft(output11_blocks(:,k));
    fft_blocks11_scale = abs(fft_blocks_double_sided11/numel(output11_blocks(:,k)));
    fft_blocks11 = fft_blocks11_scale(1:numel(output11_blocks(:,k))/2+1);
    fft_blocks11(2:end-1) = 2*fft_blocks11(2:end-1);
    average_fft_blocks1(1,k) = mean(fft_blocks11);
end

% Define block size
block_size2 = 2^10; % 4096 samples
% Calculate the number of blocks
num_blocks2 = floor(numel(recData(:,2)) / block_size2);
% Reshape the signal into a matrix with each column representing a block
output12_blocks = reshape(output1(1:num_blocks2*block_size2), block_size2, num_blocks2);
Freq_bins2 = [1:1:num_blocks2];
for k = 1:num_blocks2
    fft_blocks_double_sided12 = fft(output12_blocks(:,k));
    fft_blocks12_scale = abs(fft_blocks_double_sided12/numel(output12_blocks(:,k)));
    fft_blocks12 = fft_blocks12_scale(1:numel(output12_blocks(:,k))/2+1);
    fft_blocks12(2:end-1) = 2*fft_blocks12(2:end-1);
    average_fft_blocks2(1,k) = mean(fft_blocks12);
end


% Define block size
block_size3 = 2^11; % 4096 samples
% Calculate the number of blocks
num_blocks3 = floor(numel(recData(:,2)) / block_size3);
% Reshape the signal into a matrix with each column representing a block
output13_blocks = reshape(output1(1:num_blocks3*block_size3), block_size3, num_blocks3);
Freq_bins3 = [1:1:num_blocks3];
for k = 1:num_blocks3
    fft_blocks_double_sided13 = fft(output13_blocks(:,k));
    fft_blocks13_scale = abs(fft_blocks_double_sided13/numel(output13_blocks(:,k)));
    fft_blocks13 = fft_blocks13_scale(1:numel(output13_blocks(:,k))/2+1);
    fft_blocks13(2:end-1) = 2*fft_blocks13(2:end-1);
    average_fft_blocks3(1,k) = mean(fft_blocks13);
end

Freq_bins1 = Freq_bins1(10:(numel(Freq_bins1)-2));
Freq_bins2 = Freq_bins2(5:(numel(Freq_bins2)-2));
Freq_bins3 = Freq_bins3(5:(numel(Freq_bins3)-2));

average_fft_blocks1 =  average_fft_blocks1(10:(numel(average_fft_blocks1)-2));
average_fft_blocks2 = average_fft_blocks2(5:(numel(average_fft_blocks2)-2));
average_fft_blocks3 = average_fft_blocks3(5:(numel(average_fft_blocks3)-2));

fig2 = figure('Position', [20 20 700 700]);
subplot(3,1,1)
plot(Freq_bins1, average_fft_blocks1);
grid on; axis tight;

subplot(3,1,2)
plot(Freq_bins2, average_fft_blocks2);
grid on; axis tight;

subplot(3,1,3)
plot(Freq_bins3, average_fft_blocks3);
grid on; axis tight;

han=axes(fig2,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel("'Number of Bins'", 'FontSize',14); ylabel("Magnitude (V)", 'FontSize',14);
%--------------------------------------------------------------------------