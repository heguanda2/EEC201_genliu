function [err] = use_codebook(codebook, testing_data)
    err = 0;
    t_num = length(testing_data(:,1));
    codebook_num = length(codebook(:,1));
    for i=1:t_num
        cur_err = [];
        for j=1:codebook_num
            cur_err = [cur_err; norm(testing_data(i,:) - codebook(j,:),2)];
        end
        [val, ind] = min(cur_err);
        err = err + val;
    end
    err = err/t_num;
end