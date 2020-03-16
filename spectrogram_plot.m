clear variables;
y1 = [];
% n = length(y);
n_window = 128;
% use spectrogram() function to calculate the energy distribution
for i=7:7
    fn = sprintf('./Data/s%d.wav',i);
    figure_name = sprintf('spectrogram for speaker %d', i);
    [y,fs] = audioread(fn);
    if length(y(1,:))==1
        raw_in = y;
    else
        raw_in = y(:,1);
    end
    y = raw_in;
    y = y./max(y); % normalization here 
    y = y - mean(y);
    % y: mxn, n is the # of audio channels in the file
    %y1 = horzcat(y1,y);
    %display(size(y));
    sound(y,fs);
    %display(fn);
    %display(fs);
    %figure()
    %spectrogram(y,n_window,round(n_window*2/3), n_window,fs);
    %title(figure_name)
    pause(2);
end

dt = 1/fs;
display(dt*256);
% use custom STFT for periodgram calculation
M_incre = floor(n_window/3); 
f_digital = (1:n_window)*2/n_window;
f_analog = ((1:n_window) - n_window/2)*fs/n_window;
dt_spectrogram = M_incre * dt;
t_spectrogram = (0:floor((length(y)-n_window)/M_incre))*dt_spectrogram;
w = 0.54 - 0.46 * cos(2*pi*(1:n_window)/n_window);
%M_incre = floor(n_window/3); 
spec_total = zeros(floor((length(y)-n_window)/M_incre), n_window);
y = y - mean(y);
for i=0:floor((length(y)-n_window)/M_incre)
    y_block = y(1+M_incre*i:M_incre*i+n_window);
    y_block_w = y_block.*w';
    spec = fftshift(fft(y_block_w));
    spec_total(i+1,:) = (abs(spec).^2)./n_window;
end

figure()
surf(t_spectrogram, f_analog.*1e-3, 20*log10(spec_total'),'EdgeColor','None');
title(figure_name)
xlabel('Time (seconds)');
ylabel('Frequency (kHz)');
zlabel('Periodgram Magnitude (dB)')
ylim([-7,7])
zlim([-150,100]);
view(10,30);


