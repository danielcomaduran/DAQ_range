function [wave fwave cfs df] = wavelets_create(sampling_rate, nr_of_wavelets,scale, le)
%WAVELET creats an array of length 2 to power2 with intervals df around the
%center frequency cf with a width defined by the mode.
%wave are the wavelets in time domaine fwave in frequency domaine.
%For EMG scale = 0.3, for 2400 Hz sampling le 2^10; nr_of_wavelets = 13.
use_ck = 1;
jmax_final = nr_of_wavelets;

% le = 2^power2;
dt = 1/sampling_rate; tmax = le*dt;
df = sampling_rate/le;

%frequency and time are required for plotting
frequency = df* (0:le/2-1); frequency2 = df*(0:le/2); time = dt * (0:le-1);

%compute center frequencies and mode for Tscharner_wavelets.
%we compute one wavelet more at the low and high end of the array and
%then delete them at the end. This way the two end wavelets are not
%distorted.
jmax = jmax_final+2;
j = 1:jmax;
cfs = 1/scale*(1.45+j-2).^1.959; %array of cf in Hz.
cfs_ch = round(cfs/df);


%make matrix with f/cfs(j) in each col indexed j and compute the wavelets
inv_cfs(j) = 1./cfs(j);
diag_inv_cfs = diag(inv_cfs);
f_mat = repmat(frequency',1,jmax);
f_mat(1,:) = eps;
frcf = f_mat*diag_inv_cfs;

% Cauchy_Wavelets
    modes = scale*cfs;
    diag_modes = diag(modes);
    log_frcf = log(frcf);
    wavelet_exponent = ((-frcf+1)+log_frcf)*diag_modes;

%limit wavelet_exponent otherwied wavelets may become NaN.
if any(any(wavelet_exponent < -300))
    dum = find(wavelet_exponent < -300);
    wavelet_exponent(dum) = -300; 
end
f_wavelets = exp(wavelet_exponent);

%normalize to plateau value.
sum_operator = ones(jmax,1);
plateau = f_wavelets * sum_operator;
mean_plateau = mean(plateau(cfs_ch(3): cfs_ch(8)));
f_wavelets = f_wavelets/mean_plateau;


if use_ck
% power normalization
% use ck correction to transform the Cauchy wavelets to power_wavelets
ck = (f_wavelets.*f_wavelets)*sum_operator;
ck = ck.^-0.5;
ck(find(ck>5)) = 5;
ck(1)=1;
ck_mat = repmat(ck,1,jmax)*1/sqrt(2);
else 
    ck_mat = 1;
end
fwave = f_wavelets.*ck_mat; 
% fwave are power normalized wavelets in frequency space.

% figure(1) % show plateau of power when summing all wavelets
% clf
% y = ((fwave * sqrt(2)) .* (fwave * sqrt(2))) * sum_operator;
% plot(y)


% reduce wavelets to jmax_final, Eliminate first and last wavelet.
fwave = fwave(:,2:end-1);
cfs = cfs(2:end-1);

%With the addition of zeros the real part of wave is the real wavelet
%and the imag part the imag wavelet.
%disp(size(fwave));
wave = ifft([fwave;zeros(le/2,jmax_final)])*2; %use of matlab fft

%prepared fwave for the output. Add the zero values for the nyquist
%frequency and sqrt(2) to make the power of sum of fwave suared equal 1.
%For speed reasons matlab ifft was used which requires these adjustments.

fwave = [fwave * sqrt(2);zeros(1,jmax_final)]; 

% test = fft_real_vvt(fwave(:,13),1);
% figure(2)
% clf
% plot (test / sqrt(le))
% hold on
% plot(real(wave(:,13)),'r');


%shift wavelet to center.
wave = [wave(le/2+1:le,:);wave(1:le/2,:)];
% To subtract mean activate next line. 
% This makes a small difference for Morlet wavelets.
wave = wave - repmat(mean(wave,1),size(wave,1),1);

%% Plots
%   Enable to plot the extreme to test whether the lowest wavelet still
%   has a decent shape and whether the wavelet in time domain fits into
%   the selected range.
% figure(2)
% clf
% plot(time,5*real(wave(:,1))','b');
% hold on
% plot(time,real(wave(:,jmax_final))','r');
% figure(3)
% plot(frequency2,fwave(:,1)','b');
% hold on
% plot(frequency2,fwave(:,jmax_final)','r');
% 
% disp('wavelets created')

end %of function 

