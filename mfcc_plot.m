clear variables;
k = 26;
window_size = 256;
fn = './Data/s1.wav';
[y,fs] = audioread(fn);
%t_test = (1:length(y)) * 1/fs;
%plot(t_test,y);
%title(fn);
%xlabel('Time (s)');
%ylabel('Amplitude');
if length(y(1,:))==1
    sound_raw = y;
else
    sound_raw = y(:,1);
end
mfcc = calc_mfcc(sound_raw, fs, k, window_size);
[Y,fs] = audioread(fn);
M_incre = floor(window_size/3); 
sound_raw_block_eng = zeros(length(mfcc(:,1)),1);
for i=0:floor((length(sound_raw)-window_size)/M_incre)
    sound_raw_block_eng(i+1) = mean(abs(sound_raw(1+M_incre*i:M_incre*i+window_size)).^2);
end
%[m, max_ind] = max(sound_raw_block_eng);
%sound_raw_block_eng(max_ind) = 0;
%[m, max_ind2] = max(sound_raw_block_eng);
x = linspace(0,1,length(mfcc(:,1))) ;
y = linspace(0,1,length(mfcc(1,:))) ;
mfcc1 = mfcc(:,2:14); 
mfcc2 = zeros(length(mfcc1(:,1)),length(mfcc1(1,:)));
for i=1:length(mfcc1(:,1))
   mfcc2(i,:) = mfcc1(i,:) - mean(mfcc1(i,:)); 
end

mfcc2 = mfcc2/(max(max(abs(mfcc2))));
% highly correlated because of the overlap filterbank and FFT window

dt_tmp = floor(256/3)*1/fs;
t_tmp = (1:length(abs(mfcc2(:,1))))*dt_tmp;
figure;
surf((1:13), t_tmp, mfcc2,'EdgeColor','None');
% title()
xlabel('MFCC coefficient');
ylabel('Time (s)');
view(-85,30);
title(fn);

[my_codebook1,err] = lbg(8, 1e-2, mfcc2, 11);
figure()
plot(mfcc2(:,1), mfcc2(:,3),'bo');
hold on;
% plot(my_codebook1(:,1), my_codebook1(:,3),'r*','MarkerSize',20);
xlabel('MFCC coefficient (1)');
ylabel('MFCC coefficient (3)');
% legend('speaker1','centroid1');
title(fn);
xlim([-2 2]);
ylim([-2 2]);
% figure()
% plot(err);
% xlabel('training iterations');
% ylabel('L2 distance to centroids');



if false
    fn = './Data/s3.wav';
    [y,fs] = audioread(fn);
    %t_test = (1:length(y)) * 1/fs;
    %plot(t_test,y);
    %title(fn);
    %xlabel('Time (s)');
    %ylabel('Amplitude');
    if length(y(1,:))==1
        sound_raw = y;
    else
        sound_raw = y(:,1);
    end
    mfcc = calc_mfcc(sound_raw, fs, k, window_size);
    [Y,fs] = audioread(fn);
    M_incre = floor(window_size/3); 
    sound_raw_block_eng = zeros(length(mfcc(:,1)),1);
    for i=0:floor((length(sound_raw)-window_size)/M_incre)
        sound_raw_block_eng(i+1) = mean(abs(sound_raw(1+M_incre*i:M_incre*i+window_size)).^2);
    end
    %[m, max_ind] = max(sound_raw_block_eng);
    %sound_raw_block_eng(max_ind) = 0;
    %[m, max_ind2] = max(sound_raw_block_eng);
    x = linspace(0,1,length(mfcc(:,1))) ;
    y = linspace(0,1,length(mfcc(1,:))) ;
    mfcc1 = mfcc(:,2:14); 
    mfcc3 = zeros(length(mfcc1(:,1)),length(mfcc1(1,:)));
    for i=1:length(mfcc1(:,1))
       mfcc3(i,:) = mfcc1(i,:) - mean(mfcc1(i,:)); 
    end
    % highly correlated because of the overlap filterbank and FFT window
    figure;
    surf(mfcc3,'EdgeColor','None');
    % title()
    xlabel('MFCC coefficient');
    ylabel('Time instance');
    view(-85,30);
    title(fn);


    [my_codebook2,err] = lbg(8, 1e-2, mfcc3, 11);






    figure()
    plot(mfcc2(:,1), mfcc2(:,3),'ro');
    hold on;
    plot(my_codebook1(:,1), my_codebook1(:,3),'r*','MarkerSize',20)
    hold on;
    plot(mfcc3(:,1), mfcc3(:,3),'bs');
    hold on;
    plot(my_codebook2(:,1), my_codebook2(:,3),'b*','MarkerSize',20)
    xlabel('MFCC coefficient (1)');
    ylabel('MFCC coefficient (3)');
    legend('speaker1','centroid1','speaker2','centroid2');
    %plot(my_codebook(:,1), my_codebook(:,3),'b+')
    %title('MFCC coefficients #1 and #5 with the centroids');

end

% ret = use_codebook(my_codebook, mfcc2);


%figure()
%for i=1:k
%   plot(hm(i,:))
%   xlim([0,257]);
%   hold on;
%end