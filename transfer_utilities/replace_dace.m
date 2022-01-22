classdef replace_dace
    % this function is a replace function when cokrging
    properties
        trgx;
        trgf;
        avgf;
        name= 'replace_dacemdl';
    end
    methods
        function obj = replace_dace(x, f)
            obj.trgx = x;
            obj.trgf = f;
            obj.avgf = mean(f, 1);
        end
        
        function [f, mu] = predict(obj, x)
            f = repmat(obj.avgf, size(x, 1), 1);
            
            % Replacing the predicted values by true values if they have been evaluated
            for j=1:size(x,1)
                tmp = repmat(x(j,:), size(obj.trgx,1), 1)-obj.trgx;  % distance to trgx
                d = sum(tmp.*tmp, 2);   % distance to trgx
                id = find(d==0);   % find whether overlap with trgx
                if ~isempty(id)
                    f(j)=obj.trgf(id);                    
                end
            end
            mu = 0;
        end
    end
end
