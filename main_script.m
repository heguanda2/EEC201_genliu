clear variables;
k_mfcc = 26;
window_size = 256;
num_code = 32; % for self augmentation test, use N=8 or N=16
stepsize = 1e-2;
err_threshold = 11;
noise_var = 0;
% data augmentation, magic number 0.05, 0.1, .2
code_books = zeros(11, num_code, 13);
in_sample_err = zeros(11,1);
display_plot = false;
for i=1:11
    iter = 0;
    fn = sprintf('./Data/s%d.wav',i);
    [y,fs] = audioread(fn);
    if length(y(1,:))==1
        raw_in = y;
    else
        raw_in = y(:,1);
    end
    mfcc = calc_mfcc(raw_in, fs, k_mfcc, window_size);
    mfcc1 = mfcc(:,2:14); 
    mfcc2_test = zeros(length(mfcc1(:,1)),length(mfcc1(1,:)));
    for j_1=1:length(mfcc1(:,1))
       mfcc2_test(j_1,:) = mfcc1(j_1,:) - mean(mfcc1(j_1,:)); 
    end
    % generate more training samples with augmentation
    mfcc_final = [];
    mfcc_final = [mfcc_final; mfcc2_test];
    while iter<50
        iter = iter + 1;
        raw_in_aug = raw_in + normrnd(0, noise_var, length(raw_in), 1); % get some perturbation
        mfcc_test = calc_mfcc(raw_in_aug, fs, k_mfcc, window_size);
        mfcc1_test = mfcc_test(:,2:14); 
        mfcc2_test = zeros(length(mfcc1_test(:,1)),length(mfcc1_test(1,:)));
        for j_1=1:length(mfcc1_test(:,1))
           mfcc2_test(j_1,:) = mfcc1_test(j_1,:) - mean(mfcc1_test(j_1,:)); 
        end
        mfcc2_test = mfcc2_test/(max(max(abs(mfcc2_test))));
        mfcc_final = [mfcc_final; mfcc2_test];
    end
    [code_books(i,:,:), buf] = lbg(num_code, stepsize, mfcc_final, err_threshold);
    in_sample_err(i) = buf(end);
    if display_plot
        figure()
        plot(mfcc_final(:,1), mfcc_final(:,3),'bo');
        hold on;
        plot(code_books(i,:,1), code_books(i,:,3),'r+','MarkerSize',20);
        hold on;
        xlim([-2 2]);
        ylim([-2 2]);
        title(fn);
    end
end

% after the training, let us generate some evaluation samples
% 1st test, add gaussian noise with 0.02 variance
test_num_max = 1;
test_num = 0;
success_num = 0;
% noise_var = 0;
result_accuacy = zeros(8,1);
for i=1:8
    fn = sprintf('./TestData/s%d.wav',i);
    for j=1:test_num_max
        test_num = test_num + 1;
        err_vec = zeros(11,1);
       %  figure_name = sprintf('spectrogram for speaker %d', i);
        [y,fs] = audioread(fn);
        if length(y(1,:))==-1
            raw_in = y;
        else
            raw_in = y(:,1);
        end
        raw_in_aug = raw_in + normrnd(0, noise_var, length(raw_in), 1); % add the perturbation here
        mfcc_test = calc_mfcc(raw_in_aug, fs, k_mfcc, window_size);
        mfcc1_test = mfcc_test(:,2:14); 
        mfcc2_test = zeros(length(mfcc1_test(:,1)),length(mfcc1_test(1,:)));
        for j_1=1:length(mfcc1_test(:,1))
           mfcc2_test(j_1,:) = mfcc1_test(j_1,:) - mean(mfcc1_test(j_1,:)); 
        end
        mfcc2_test = mfcc2_test/(max(max(abs(mfcc2_test))));
        if display_plot
            figure()
            plot(mfcc2_test(:,1), mfcc2_test(:,3),'bo');
            title(fn);
            xlim([-2 2]);
            ylim([-2 2]);
        end
        % hold on;
        for k=1:11
            err_vec(k) = use_codebook(squeeze(code_books(k,:,:)), mfcc2_test);
        end
        [val, ind] = min(err_vec);
        if ind==i
            success_num = success_num + 1;
        end
    end
    result_accuacy(i) = success_num/test_num;
    test_num = 0;
    success_num = 0;
    display(fn);
    display(err_vec');
end
display(result_accuacy');
