function [mfcc,sound_raw] = calc_mfcc(sound_raw,fs, k,window_size)
    % fn: directory and file name of the .wav file 
    % k: filter bank order 
    % window_size: FFT window size
    % fs = 12500;
    % [sound_raw,fs] = audioread(fn);
    % normalization
    sound_raw = sound_raw ./ max(sound_raw);
    % for test purpose
    % fs = 12.5e+3;
    % sound_raw = sin((1:13056).*(1/fs)*pi*500)';
    % dt = 1/fs;
    M_incre = floor(window_size/3); 
    % f_digital = (1:window_size)*2/window_size; % unit:pi
    % f_analog = ((1:window_size) - window_size/2)*fs/window_size; %unit:Hz
    hamming_w = 0.54 - 0.46 * cos(2*pi*(1:window_size)/window_size);
    % spec_total = zeros(floor((length(sound_raw)-window_size)/M_incre), window_size);
    spec_total = zeros([] , window_size);
    ind_use = 0;
    sound_raw_block_true = [];
    for i=0:floor((length(sound_raw)-window_size)/M_incre)
        sound_raw_block = sound_raw(1+M_incre*i:M_incre*i+window_size);
        sound_raw_block_eng = mean(abs(sound_raw_block).^2);
        if sound_raw_block_eng > 1e-4
            sound_raw_block_w = sound_raw_block.*hamming_w';    
            spec = fftshift(fft(sound_raw_block_w)); 
            spec_total = [spec_total; (abs(spec').^2)/window_size]; % periodgram = abs(fft)^2
            ind_use = ind_use + 1;
            sound_raw_block_true = [sound_raw_block_true; sound_raw_block_w];
        end
        % sound_raw_block_w = sound_raw_block.*hamming_w';
        % spec = fftshift(fft(sound_raw_block_w)); 
        % spec_total(i+1,:) = (abs(spec).^2)/window_size; % periodgram = abs(fft)^2
    end
    % figure()
    % plot(sound_raw_block_true);
    
    dt_tmp = M_incre*1/fs;
    t_tmp = (1:length(abs(spec_total(:,1))))*dt_tmp;
    df_tmp = 12.5/256; 
    f_tmp = (-128:127)*df_tmp;
    %figure()
    %surf(f_tmp, t_tmp, 20*log10(abs(spec_total)),'EdgeColor','None');
    %title('periodgram before mel FB for s1');
    
    fb = melfb(k, window_size, fs);
    %figure()
    for m1 = 1:k
        if m1>k-10
            fb(m1,:) = fb(m1,:) * 1.3;
        end
        %plot(f_tmp,fb(m1,:));
        %hold on;
    end
    % spec_flt = zeros(length(spec_total(:,1)), window_size); % mel-scale spectrum 
    spec_flt_histo = zeros(length(spec_total(:,1)), k); % mel-scale spectrum 
    for i=1:length(spec_total(:,1))
       % buf = zeros(1,window_size); 
       for j = 1:k
          % buf = buf + fb(j,:).*spec_total(i,:);  
          spec_flt_histo(i,j) = dot(fb(j,:), spec_total(i,:));
       end
       % spec_flt(i,:) = buf; 
    end
    
    % figure()
    % surf((1:26), t_tmp, 20*log10(abs(spec_flt_histo)),'EdgeColor','None');
    % title('periodgram after mel FB');
    % do a test here
    % i = 90;
    % figure()
    % plot(spec_total(i,:));
    % figure()
    % plot(spec_flt_histo(i,:));
    % we should notice spec_flt_histo is a smoothed version of spec_total
    
    
    % now go to cepstrum
    % spec_flt_log = log10(spec_flt);
    spec_flt_log = log10(spec_flt_histo);
    
    mfcc = my_dct(spec_flt_log,0);
    
    %sound(y,fs);
    %pause(5);
end
