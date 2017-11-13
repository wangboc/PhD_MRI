function ImgRecon = SENSEArbitrary_regul_GN(kspace_data, sen_map, InitImg, trajectory,...
                                        Density, Ls, afa, inter_num)

%transform the matrix KspaceDataWeighted which has the size nc*nk1 to vector m 
%nk1 is less than nk which is original number of K-space data because KspaceDataWeighted is down-sampled data 

%calculate the diagnoal matrix I which has the size N^2*N^2
[D1, D2, Coilnum] = size(sen_map);

weight = sqrt(sum(abs(sen_map) .^ 2, 3));
I=ones(D1, D2) ./ (weight + eps);

%calculate the diagnoal matrix D which has the size (nc*nk1)*(nc*nk1)
%KspaceDataWeighted=FullKspaceData(:,(off+1):ReduceFactor:(Dk-off));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CG algorithm
%RequiredAcc=1/100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%a=I*E(H)*D*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s = 1 : Coilnum
    xspaceimage(:, :, s) = SEH_GN(trajectory, kspace_data(:, :, s), Density);
    xspaceimageWeighted(:, :, s) = xspaceimage(:, :, s) .* conj(sen_map(:, :, s));
end
%load my_data3 %just load xspaceimageWeighted

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = I .* sum(xspaceimageWeighted, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 以下采用共轭梯度算法可在SENSE论文附件C中查看
% Pruessmann K P, Weiger M, B?rnert P, et al. 
% Advances in sensitivity encoding with arbitrary k‐space trajectories[J]. 
% Magnetic resonance in medicine, 2001, 46(4): 638-651.
RequiredAcc = 1e-10;
b = InitImg;
p = a;
r = a;
count = 0;
delta = r(:)' * r(:) / (a(:)' * a(:));

while (count < inter_num && delta > RequiredAcc)
    
    count = count + 1;
        q = SMatrixAfunction3_GN(p, sen_map, trajectory, I, Density); 
        Lq = Ls * reshape(p, D1 * D2, 1);
      %  Lq=Ls*reshape(q,n1*n2,1);
%         bs = sum(sum(abs(b)));
%         if bs ~= 0
%           % bs=bs/((n1*n2)*1.18^count);
%            bs = bs / (D1 * D2);
%         else
%            bs = 1;
%         end
        %q=q+afa*reshape(Lq,n1,n2)/bs;
        q = q + afa * reshape(Lq, D1, D2);
        %call function MatrixAfunction to calculate  I*E(H)*D*E*I
        b = b + (r(:)' * r(:)) / (p(:)' * q(:)) * p;
        
        tempr = r(:);
        r = r - (r(:)' * r(:)) / (p(:)' * q(:)) * q;
        p = r + (r(:)' * r(:)) / (tempr' * tempr) * p;
     %   p=I.*p; %%%
        delta = r(:)' * r(:) / (a(:)' * a(:));
        %figure, imshow(abs(rot90(I.*b,1)),[])
end
%ImgRecon=b;        
ImgRecon = I .* b;







