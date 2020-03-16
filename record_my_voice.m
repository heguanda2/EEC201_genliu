clear variables;

recObj = audiorecorder(12500, 8, 1);
disp('Start speaking.')
recordblocking(recObj, 5);
disp('End of Recording.');

y = getaudiodata(recObj);

plot(y);
% audiowrite('s12.wav', y, 12500);

% audiowrite('s12_test.wav', y, 12500);