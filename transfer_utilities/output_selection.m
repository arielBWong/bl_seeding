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
    
    front1_X(corner_indx, :) = [];                            % X needs no normalization
    %---------------------------------------
    
    % Normalization process for the third point
    [corners_min, ~] = min(front1_F, [], 1);
    [corners_max, ~] = max(front1_F, [], 1);
    
    
    % normalize
    front1_Fnorm = (front1_F - corners_min) ./ (corners_max - corners_min);
    corner1_norm = front1_Fnorm(corner_indx(1), :);
    corner2_norm = front1_Fnorm(corner_indx(2), :);
    corners_norm = [corner1_norm; corner2_norm];
    
    % remove corners in F
    front1_Fnorm(corner_indx, :) = [];
    rest_frontnorm = front1_Fnorm;
    
    
    
    
    % find middle point of normalized
    d = pdist2(rest_frontnorm, corners_norm,  'euclidean');
    d = abs(d(:, 1) - d(:, 2));   % only two corners
    [~, ord] = sort(d);
    third_idx = ord(1);
    
    %     % seperate the rest
    %     front1_F(corner_indx, :) = [];
    %     rest_front = front1_F;
    %     front1_X(corner_indx, :) = [];
    %
    %     % find the middle
    %     d = pdist2(rest_front, corners, 'euclidean');
    %     d = abs(d(:, 1) - d(:, 2));   % only two corners
    %     [~, ord] = sort(d);
    %     third_idx = ord(1);
    
    %----------------------
    % collect output
    middle_norm = rest_frontnorm(third_idx, :) ;
    middle = corners_min + middle_norm .* (corners_max - corners_min);   % denormalization
    best_f = [corners;  middle];
    best_x = [cornersx; front1_X(third_idx, :)];  % front1_X checked removed cornerX
    
    
else
    % only choose corners, because only the first front went through crowding
    % sort
    
    front1_F = pop.F;
    front1_X = pop.X;
    
    % pick up two corners
    [~, corner_indx] =  min(front1_F, [], 1);
    corner1 =  front1_F(corner_indx(1), :);
    corner2 =  front1_F(corner_indx(2), :);
    
    corner1x =  front1_X(corner_indx(1), :);
    corner2x =  front1_X(corner_indx(2), :);
    
    
    best_f = [corner1;  corner2];
    best_x = [corner1x; corner2x];
end

end
