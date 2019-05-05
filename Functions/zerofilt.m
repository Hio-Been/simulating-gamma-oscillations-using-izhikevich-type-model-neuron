%ZEROFILT zero-phase bandpass filtering with firls.m and filtfilt.m
%   result = zerofilt(data,hp_freq,lp_freq,f_samp);
%   For HP filtering, fill lpfreq with Nyq. 
%   For LP filtering, fill hpfreq with 0.
%   Data length must be bigger than 1200.
%   Sampling frequency fs is recommended to be smaller than 200 Hz.
%   Otherwise filtering performance decreases.

function result = zerofilt(data,hpfreq,lpfreq,fs)

Nyq = fs/2;

if hpfreq == 0,
    % low-pass (when hpfreq = 0)
    f = [ 0 lpfreq lpfreq Nyq]/Nyq;
    a= [ 1 1 0 0 ];
    b = firls(500,f,a);     % least-square linear phase FIR filter of order 400
elseif lpfreq == Nyq,
    % high-pass (when lpfreq = Nyq);
    f = [ 0 hpfreq hpfreq Nyq]/Nyq;
    a= [ 0 0 1 1];
    b = firls(500,f,a);     % least-square linear phase FIR filter of order 400
else 
    f = [ 0 hpfreq hpfreq lpfreq lpfreq Nyq]/Nyq;
    a= [ 0 0 1 1 0 0];
    b = firls(500,f,a);     % least-square linear phase FIR filter of order 400
end

% test filter design
% [H,f] = freqz(b,1,512,2);
% figure; plot(f,abs(H));

result = filtfilt(b,1,data);        % zero-phase filtering







