function [hm] = melfb(p,n,fs)
    % p: # of filters in filterbank
    % n: length of FFT
    % fs: sampling rate
    f_low = 50;
    f_high = fs/2;
    f_mel_low = 1125 * log(1+f_low/700);
    f_mel_high = 1125 * log(1+f_high/700);
    f_mel_freq_ptrs =(0:p+1)*(f_mel_high - f_mel_low)/(p+1) + f_mel_low;
    f_fb_ptrs = 700 * (exp(f_mel_freq_ptrs/1125) - 1);
    f_bin_round = floor(n*f_fb_ptrs/fs);
    hm = zeros(p, n);
    for filter_num=1:p
        for k=1:n/2
            if k<f_bin_round(filter_num)
                hm(filter_num, k+n/2)= 2e-22; % very small value to avoid -inf after log
            elseif f_bin_round(filter_num) <= k && k < f_bin_round(filter_num+1)
                hm(filter_num,k+n/2) = (k-f_bin_round(filter_num))/(f_bin_round(filter_num+1)-f_bin_round(filter_num));
            elseif f_bin_round(filter_num+1) <= k && k < f_bin_round(filter_num+2)
                hm(filter_num,k+n/2) = (f_bin_round(filter_num+2)-k)/(f_bin_round(filter_num+2)-f_bin_round(filter_num+1));            
            else
                hm(filter_num, k+n/2)= 2e-22; % very small value to avoid -inf after log
            end
        end
        hm(filter_num,1:n/2)=flip(hm(filter_num,n/2+1:end));
    end
end