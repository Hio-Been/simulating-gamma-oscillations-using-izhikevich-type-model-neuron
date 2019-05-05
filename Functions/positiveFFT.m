function [X,freq] = positiveFFT(x,Fs,plotOption)
%%  [X,freq]=positiveFFT(x,Fs,plotOption)

N=(length(x));  %get the number of points
k=0:N-1;        %create a vector from 0 to N-1
T=N/Fs;         %get the frequency interval
freq=k/T;       %create the frequency range
cutOff = ceil(N/2);
freq = freq(1:cutOff);


nTrials = size(x, 1);
if nTrials == 1
    X=fft(x)/N*2; % normalize the data
%     disp(size(X))
    %only want the first half of the FFT, since it is redundant
    X = X(1:cutOff); % Single trial FFT
else
    X = [];
    for trialIdx = 1:nTrials
        X_temp = fft(x(trialIdx,:))/N*2;
        size(X_temp)
        X_temp = X_temp(1:cutOff); % Single trial FFT
        X = [X; X_temp];
    end
end

if plotOption
    for trialIdx = 1:nTrials
        plot(freq, abs((X(trialIdx,:))));
        hold on;
    end;
    xlim([1,161])
    xlabel('Freq (Hz)')
    ylabel('Normalized Amplitude')
    
    hold off;
end

return
