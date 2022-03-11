function [best_x, best_f, best_c] = output_selection(fronts, pop)
% This function is used  for gsolver to return solutions 
best_c = [];
nf1 = length(fronts(1).f);

if nf1 >= 3
    indx1 = fronts(1).f;
    front1_F = pop.F(indx1, :);
    front1_X = pop.X(indx1, :);
    
    % pick up two corners
    [~, corner_indx] =  min(front1_F, [], 1);
    corner1 =  front1_F(corner_indx(1), :);
    corner2 =  front1_F(corner_indx(2), :);
   
    corner1x =  front1_X(corner_indx(1), :);
    corner2x =  front1_X(corner_indx(2), :);
    
    corners = [corner1; corner2];
    cornersx=  [corner1x; corner2x];
    
    % seperate the rest
    front1_F(corner_indx, :) = [];
    rest_front = front1_F;
    front1_X(corner_indx, :) = [];
    
    % find the middle
    d = pdist2(rest_front, corners, 'euclidean');
    d = abs(d(:, 1) - d(:, 2));   % only two corners
    [~, ord] = sort(d);
    third_idx = ord(1);
    

    best_f = [corners;  rest_front(third_idx, :)];
    best_x = [cornersx; front1_X(third_idx, :)];
    

else
    % only choose corners, because only the first front went through crowding
    % sort 
    best_f = pop.F(fronts(1).f, :);
    best_x = pop.X(fronts(1).f, :);    
end

end
