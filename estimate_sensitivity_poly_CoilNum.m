function  WeightingFunctions=estimate_sensitivity_poly_CoilNum(row, column, CoilNum, x0, y0, s_poly, order)

%  ***********************************************
%  Written by Jinhua Sheng at 09/27/2007 
%  University of Wisconsin
%  x0 for averaged value of x
%  y0 for averaged value of y
%  s_poly contains the coefficients calculated by polynomial fitting
%  **********************************************
for s = 1 : CoilNum
    A_poly(1:order+1, 1:order+1) = s_poly(1:order+1, 1:order+1, s);
    for ku = 1 : row 
        for kv = 1 : column
            P(1:order+1, 1) = (ku / row - x0) .^ (0:order);
            Q(1, 1:order+1) = (kv / column - y0) .^ (0:order);
            PQ = P * Q;
            SPQ = A_poly .* PQ;
            WeightingFunctions(ku, kv, s) = sum(sum(SPQ(:, :)));
        end
    end
  %  WeightingFunctions(:,:,s)=WeightingFunctions(:,:,s);%.*WFT;
end

