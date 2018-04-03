%% DAQ_EMG_range                                                02/04/2018
%  - This script takes the data from the biofeedback experiment [Comaduran
%    2016] and down-samples the data dynamic range to the values in 
%    src.dynamic
%  - The output of the script are:
%    - A graph that shows how the raw coherence changes with the changing
%      dynamic range
%    - The 'h' and 'p' values for a paired t-test comparing the coherence
%      of interest between control and biofeedback conditions

%% Import data
load('raw_EMG_event');

%% Variables
src.rate = 2400;                    % Sampling frequency [Hz]
src.subjects = 1:10;                % Number of subjects
src.trials = 1:6;                   % Number of trials
src.bin = 2^11;                     % Number of point for FFT calculation and peak finding
src.amp;                            % Matrix with amplification factor 
src.range = [-2, 2];                % Original dynamic range of the sampled data
src.dynamic = [14, 12, 10, 8, 6, 4, 2];      % Resolution in [Bits]
% src.dynamic = [14, 1];      % Resolution in [Bits]
                                    % - (Original, downsamples...)
                                    % - SENIAM recommended [12Bit]

w = rectwin(src.bin);   % Number of  data points to obtain X segments
n = 0;                  % Number of points to overlap
L = 90;                 % Number of segments for analysis 
a = 0.95;               % Confidence interval for Rosenberg significance

% Preallocate variables
Cxy = cell(length(src.subjects), length(src.trials));       % Concurrent coherence [Frequency, Coherence]
CoI = zeros(length(src.subjects), length(src.trials));      % Concurrent coherence of interest [%]
h = zeros(length(src.dynamic), 1);
p = zeros(length(src.dynamic), 1);

% Rosenberg significance for coherence calculation
ros = 1 - (1-a)^(1/(L-1));

%% Data analysis
% Frequency vector - for Coherence and PSD
T = src.bin * (1/src.rate); % Period
df = 1 / T;                 % Frequency differential
freq = (0:src.bin/2)' * df; % Frequency vector

figure
hold on;

for d = 1:length(src.dynamic)
    for s = src.subjects    
        for t = src.trials
            src.d = d;
                        
            %% Change dynamic resolution
            measured_data = data_event{s,t}*src.amp(s,t);
            resampled_data = down_range(src, measured_data); 
            
            %% Calculate coherence
            Cxy{s,t} = mscohere(resampled_data(:,1), resampled_data(:,2), w, n, src.bin, src.rate);   % Simultaneous events
            
            %% Normalize to coherence of interest
            fi = [find(freq > 10,1,'first'), find(freq < 100,1,'last')];
            CoI(s,t) = trapz(freq(fi(1):fi(2)),  Cxy{s,t}(fi(1):fi(2))) / trapz(freq(fi(1):fi(2)), ones(size(Cxy{s,t}(fi(1):fi(2)))));  % Concurrent
            
        %% Calculate intermuscular coherence
%         z = circshift(data_event(:,2), src.bin);    % Shift data by src.bin points
            
%         Cxz{s,t} = mscohere(data_event(:,1), z, w, n, src.bin, src.rate);                 % Consecutive events      
%         
%         %% Calculate intermuscular coherence with randomized VM
%         m = length(data_event)/src.bin;         % # of columns to reorder
%         i = randperm(m);                        % Randomize columns
%         c = reshape(data_event(:,2), [], m);    % Reshape vector to matrix
%         r = reshape(c(:,i), [], 1);             % Reshape matrix to vector
%         Rxy{s,t} = mscohere(data_event(:,1), r, w, n, src.bin, src.rate);   % Randomized coherence
%                 

%         SoI(s,t) = trapz(freq(fi(1):fi(2)),  Cxz{s,t}(fi(1):fi(2))) / trapz(freq(fi(1):fi(2)), ones(size(Cxz{s,t}(fi(1):fi(2)))));  % Subsequent
%         RoI(s,t) = trapz(freq(fi(1):fi(2)),  Rxy{s,t}(fi(1):fi(2))) / trapz(freq(fi(1):fi(2)), ones(size(Rxy{s,t}(fi(1):fi(2)))));  % Random
%                         
        end
    end
    
    %% Plot coherence
    all_Cxy = cell2mat(reshape(Cxy, 1, []));
    hold on;
    plot(freq, mean(all_Cxy, 2), 'DisplayName', [num2str(src.dynamic(d)), ' bits'])

    %% T-test = 
    x = reshape(CoI(:,1:3), [], 1);
    y = reshape(CoI(:,4:6), [], 1);
    [h(d), p(d)] = ttest(x,y);
    
%     disp(' '); disp(['Subject ' num2str(s)]);
%     disp(['Trial ' num2str(t) ' - ' num2str(n_events) ' events found']);
end

hold off;
toc

%% Extra functions
function output = down_range(src, data)
    %% down_range
    %  This function takes each column in 'data', digitizes it by n-bits as
    %  stated in src.dynamic(1), and then returnts the resampled data by
    %  src.dynamic(2:end)
    
    output = zeros(size(data));
     
    for i = 1:size(data,2)
        analog_original = reshape(data(:,i), src.bin, []) + src.range(2);
        digital_original = floor((analog_original*(2^src.dynamic(1))/diff(src.range)));
        
        digital_new = floor(digital_original*(2^src.dynamic(src.d))/(2^src.dynamic(1)));
        analog_new = digital_new*diff(src.range)/(2^src.dynamic(src.d)) - src.range(2);
        
        output(:,i) = reshape(analog_new, [], 1);
    end
    
%     if src.d == 2
%         figure
%         hold on;
%         plot(analog_original(:,1))
%         plot(analog_new(:,1))
%     end
end

function [output, n_events, pow_emg, time_emg] = eventDetection(src, data)
    %% eventDetection
    %  - This function takes a column vector of EMG activity and outputs
    %  all thelocations of EMG peaks found for further coherence analysis.
    % - data(:,1) = no use
    % - data(:,2) = used for peak finding
    % - output = concatenated data from each peak +/- (src.bin/2)
    
    %% Append zeros to speed-up wavelet calculation
    trimPow2 = nextpow2(length(data));                  % Calculate next power of 2
    zeroPad = 2^trimPow2 - length(data);                % Calculate zeros to add
    data_zeropad = padarray(data, zeroPad, 0, 'post');  % Add zeros to data
    
    %% Perform wavelet transform
    le = length(data_zeropad);                                                  % Length of mother wavelet
    nw = 11;                                                                    % Number of wavelets
    [~, fftWave, ~, ~] = wavelets_create(src.rate, nw, 0.3, le);                % Create wavelets
    
    waveletConv = fftWave .* repmat(fft_real_vvt(data_zeropad(:,2)', 0), 1, nw);% Wavelet transform for VM
    ifftWaveletConv = fft_real_matrix(waveletConv, 1);                          % - Calculate FFT of wavelets
    power = (ifftWaveletConv .* conj(ifftWaveletConv))';                        % - Autocorrelation for power
    waveletsNo = 3:11;                                                          % - Select wavelets for envelope
    envelope_vm = sum(sqrt(power(waveletsNo, :)), 1);                              % - Obtain envelope
    
        
    wavelet_vl = fftWave .* repmat(fft_real_vvt(data_zeropad(:,1)', 0), 1, nw);% Total power per event VM
    ifftWaveletConv = fft_real_matrix(wavelet_vl, 1);                          
    power_vl = (ifftWaveletConv .* conj(ifftWaveletConv))';
    envelope_vl = sum(sqrt(power_vl(waveletsNo, :)), 1);
    
    
    % Reshape envelope to interval bins
    interval = 256;                                 
    envelopeRS = reshape(envelope_vm, interval, []);
    envelope2 = mean(envelopeRS, 1);
    envelope2 = gliding_filter(envelope2, 5);
    
    %% Find midpoint for the event
%     tAmp = 0.45; % [%] of Max EMG activation to find peak
    [~, LOCS] = findpeaks(envelope2, 'MinPeakHeight', max(envelope2) * src.amp);
    peak_locations = LOCS * interval - interval / 2;
    
    %% Output
    % Event data
    n_events = numel(peak_locations) - 2;
    bit_mat = ones(length(peak_locations)-2, src.bin) .* ((-src.bin/2):1:(src.bin/2)-1);    % Matrix to calculate active EMG bins
    output_idx = reshape((bit_mat + ones(size(bit_mat)).*peak_locations(2:end-1)')',[],1);   % Vector of indices used for output
    
    output = zeros(length(output_idx), 2);  % Use output_idx for selecting data points
    output(:,1) = data(output_idx, 1);      % - VL
    output(:,2) = data(output_idx, 2);      % - VM 
    
    % Power data by EMG event
    pow_emg(:,1) = sum(reshape(envelope_vm(output_idx), src.bin, []), 1)';  % Total power per event [VL]
    pow_emg(:,2) = sum(reshape(envelope_vl(output_idx), src.bin, []), 1)';  % Total power per event [VM]
    
    % Time between events
    time_emg = diff(peak_locations/src.rate)';
    
end

