function [code_book, err_final] = lbg(M_req, stepsize,train_data, err_threshold)
    % we should have one code book for each user
    % train_data TxK matrix, T is the time instance, K is the 13 MFCC
    % coefficients
    code_book = [];
    tmp_book_old = [];
    tmp_book_new = [];
    M = 1;
    cur_centroid = mean(train_data); 
    cur_err = 0;
    t_num = length(train_data(:,1));
    feature_num = length(train_data(1,:));
    err_final = [];
    for i=1:t_num
        cur_err = cur_err + norm(train_data(i,:)-cur_centroid, 2);
    end
    cur_err = cur_err/t_num;
    if cur_err < err_threshold && M_req==1
       code_book = [code_book; cur_centroid]; 
    else
        tmp_book_old = [tmp_book_old; cur_centroid];
        while M<M_req
            M = 2*M;
            tmp_book_new = [];
            for j=1:M/2
                buf_new_1 = tmp_book_old(j, :) * (1 + stepsize);
                buf_new_2 = tmp_book_old(j, :) * (1 - stepsize);
                tmp_book_new = [tmp_book_new; buf_new_1];
                tmp_book_new = [tmp_book_new; buf_new_2];
            end
            % finish the splitting
            tmp_book_old = tmp_book_new; 
            cur_err = 1000; % assign a high value to init the loop
            while cur_err > err_threshold
               % M = 2*M;
               new_cur_err = 0;
               new_centroid = zeros(M, feature_num);   
               M_sample_cnt = zeros(1,M); % # of training vector assigned to each cell
               for i=1:t_num
                   err_to_each_code = [];
                   for j=1:M
                       err_to_each_code = [err_to_each_code; norm(train_data(i,:)-tmp_book_old(j,:), 2)];
                   end
                   [val, ind] = min(err_to_each_code);
                   new_centroid(ind,:) = new_centroid(ind,:) + train_data(i,:);
                   M_sample_cnt(ind) = M_sample_cnt(ind) + 1;
                   new_cur_err = new_cur_err + val;
               end
               for i=1:M
                  if(M_sample_cnt(i)>0)
                     new_centroid(i,:) = new_centroid(i,:)/(M_sample_cnt(i)); 
                  end
               end
               % update the codebook
               tmp_book_old = new_centroid; 
               new_cur_err = new_cur_err/t_num;
               cur_err = new_cur_err;
               err_final = [err_final; cur_err];
               display(err_final);
            end
        end
        code_book = [code_book; tmp_book_old];
    end
    % err_final = cur_err;
end