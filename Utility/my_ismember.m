 function member_idx = my_ismember(A, B)
 % can deal with A has repeated points and B also have repeated point
 n                    = size(A, 1);
 member_idx           = [];
 for i = 1:n
    A_row             = A(i, :);
    [in_B, row_idx]   = ismember(A_row, B, 'rows');
    if in_B
        B(row_idx, :) = B(row_idx, :) + 1;
    end
    member_idx        = [member_idx, row_idx];
 end