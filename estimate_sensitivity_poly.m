function  WeightingFunction=estimate_sensitivity_poly(row,column,A_poly,order)

%  ***********************************************
%  Written by Jinhua Sheng at 09/27/2007 
%  University of Wisconsin
%  **********************************************
x0=sum(1:row)/row^2;
y0=sum(1:column)/column^2;

for ku = 1 : row 
    for kv = 1 : column
        P(1:order+1,1)=(ku/row-x0).^(0:order);
        Q(1,1:order+1)=(kv/column-y0).^(0:order);
        PQ=P*Q;
        SPQ=A_poly.*PQ;
        WeightingFunction(ku,kv)=sum(sum(SPQ(:,:)));
    end
end

