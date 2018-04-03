%% This function filters the 60 noise as well as the specified BPF
%   - x1 is entered as a column vector
%   - For the time being it is expected that the signal is 2^13 long this is about 3.41 sec

function [filtered_sig] = EMG_line_filter(sig, power_of_2, line_frequency, test)
%% Global variables
N = 2^power_of_2;
N4 = 2^7;
samplingrate = 2400;
time = ((1:N)-1)*1/samplingrate;
nyqist_frequency = N/2+1;
T = N/samplingrate;
df = 1/T;

[filter1, filter60Hz, frequency_range] = EMG_filter(power_of_2, samplingrate, line_frequency, test);

%% Iterate through the signal
% Reorient data into column vector
if size(sig,1)== 1
    sig = sig';
end

filtered_sig = sig;
sig_length = length(sig);
flag_on = true;
flag_last_signal = false;
pointer = 1;

while flag_on
    y = sig(pointer:pointer-1 + N);
    
    %% Apply filters
    fftx = fft_real_vvt(y,0); %returns a col_vector
    fftx2 = fftx' .* filter1; % remove 400 Hz and very low frequency
    x2 = fft_real_vvt(fftx2,1)';
    fftx_filtered = fftx' .* filter60Hz; %xFiltered has only the signal from 60 Hz and its harmonics
    xFiltered = fft_real_vvt(fftx_filtered,1)'; %is a row vector
    
    %% Partition filter
    intervall_length = round((1/line_frequency)* samplingrate);
    intervalls = floor(N/intervall_length);
    Smatrix = reshape(xFiltered(1:intervall_length*intervalls), intervall_length,intervalls); %This is shorter than whole file
    Smatrix2 = reshape(x2(1:intervall_length*intervalls),intervall_length,intervalls);
    
    m = mean(Smatrix,2);
    m = m/norm(m);
    
    if test
        figure(16)
        clf
        plot(m)
    end
    
    M = repmat(m,1,intervalls+1);
    
    weights = m' * Smatrix2;
    meanweights = mean(weights);
    weights(end+1) = 0;
    weights = weights * 0 + meanweights;
    weights = diag(weights);

    S60HZ = M*weights;
    S60HZ = reshape(S60HZ,1,intervall_length *(intervalls+1));
    S60HZ = S60HZ(1:N);
    
    Scleaned = x2' - S60HZ';
    filtered_sig(pointer:pointer-1 + N) = Scleaned - mean(Scleaned);
    
    if flag_last_signal
        flag_on = false;
    else
        if pointer-1 + N + N >= sig_length
            flag_last_signal = true;
            pointer = sig_length - N + 1;
        else
            pointer = pointer + N;
        end
    end
end
end

function [filter1 filter60Hz frequency_range] = EMG_filter(power_of_2, samplingrate, line_frequency, test)
%Creats the EMG filter
%   2 to the power_of_2 indicates length of the accepted signal, recomended
%   12 or 13. 

% Define basic variables
N = 2^power_of_2;
time = ((1:N)-1)*1/samplingrate;
nyqist_frequency = N/2+1;
T = N/samplingrate;
df = 1/T;

% Define frequency range
FP = (1:nyqist_frequency);
FP2 = (1:nyqist_frequency-1)*df; % leave out the first point
frequency_range = (FP-1)*df; % Frequency in Hertz

% line frequency extractor filter
% these are band filters for multiples of 60 Hz
CF1 = line_frequency; % line frequency
factor = 2;
F1 = FP2/CF1; F2 = FP2/(CF1*2);
F3 = FP2/(CF1*3);F4 = FP2/(CF1*4);F5 = FP2/(CF1*5);
filter = exp((-F1+1+log(F1))*factor*CF1)...
    + exp((-F2+1+log(F2))*factor*(CF1*2))...
    + exp((-F3+1+log(F3))*factor*(CF1*3))...
    + exp((-F4+1+log(F4))*factor*(CF1*4))...
    + exp((-F5+1+log(F5))*factor*(CF1*5));
filter60Hz = [0 filter];

% Low pass filter 500 Hz
flowpass = 500;
factor = 0.3;
F6 = FP2/(flowpass);
filter6 = exp((-F6+1+log(F6))*factor*(flowpass));
filter6 = [0 filter6];
%filter6 = 1-filter6; %take nothing above flowpass
limit = floor(flowpass/df);
filter6(limit:end) = 1;

% High pass filter, for EMG 15 Hz, second wavelet
f_high_pass = 15;
factor = 0.3;
F7 = FP2/(f_high_pass);
filter7 = exp((-F7+1+log(F7))*factor*(f_high_pass));
filter7 = 1-[0 filter7];%Take everything above 18 Hz
limit = floor(f_high_pass/df);
filter7(limit:end) = 0;

% Band supression 400 Hz
f400 = 400;
factor = 20;
F8 = FP2/(f400);
filter8 = exp((-F8+1+log(F8))*factor*(f400));
filter8 = [0 filter8] * 0.99; %0.99 to prevent from crshing when computing coherence

filter1 = 1 - (filter6 + filter7 + filter8);

if test == 1
    figure(10)
    clf
    plot(frequency_range,filter60Hz, 'g')
    hold on
    plot(frequency_range,filter6,'b')
    plot(frequency_range,filter7+0.1,'r')
    plot(frequency_range,filter8,'g')
    plot(frequency_range,filter1,'k','LineWidth',3)
end
end