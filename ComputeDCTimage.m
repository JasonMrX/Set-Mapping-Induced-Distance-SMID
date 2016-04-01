function C = ComputeDCTimage(X, T, type)
% type = 1 : forward
%      = 2 : inverse
    [row, col] = size(X) ;
    N = size(T, 1);
    if type == 2
        T = T';
    end
    for j = 1 : N : row
        for i = 1 : N : col
            block = X(j : j + N - 1, i : i + N - 1) ;
            C(j : j + N - 1, i : i + N - 1) = (T * block) * (T') ;
        end
    end
end