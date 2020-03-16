function [out] = my_dct(in, dim)
    % in: nxm matrix
    % dim: determine wheter we compute dct along column or row
    % when dim = 0, compute along row
    % when dim = 1, compute along column
    row_len = length(in(1,:));
    col_len = length(in(:,1));
    % we might consider doing an exception catch here in case in is a
    % vector
    out = zeros(col_len, row_len);
    if dim==0
        for i=1:col_len
           buf = in(i,:);
           for j=1:row_len
              vec_cos = cos((j-1)*((1:row_len) - 0.5)*pi/row_len);
              out(i,j) = dot(buf, vec_cos); 
           end
        end
    else
        for i=1:row_len
           buf = in(:,i);
           for j=1:col_len
              vec_cos = cos((j-1)*((1:col_len) - 0.5)*pi/col_len);
              out(j,i) = dot(buf, vec_cos); 
           end
        end
    end
end