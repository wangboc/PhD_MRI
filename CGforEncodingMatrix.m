function [ WeightingFunctions ] = CGforEncodingMatrix(kspace_data, Image, sen_map, InitImg, trajectory,...
                                        Density, Ls, afa, inter_num )
%CGFORENCODINGMATRIX 此处显示有关此函数的摘要
%   此处显示详细说明


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
    xspaceimage(:, :, s) = ifft2(fftshift(kspace_data(:, :, s)));
    xspaceimageWeighted(:, :, s) =  conj(Image) .* (xspaceimage(:, :, s));
end
%load my_data3 %just load xspaceimageWeighted

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 以下采用共轭梯度算法可在SENSE论文附件C中查看
% Pruessmann K P, Weiger M, B?rnert P, et al. 
% Advances in sensitivity encoding with arbitrary k‐space trajectories[J]. 
% Magnetic resonance in medicine, 2001, 46(4): 638-651.
for s = 1 : Coilnum
    a = xspaceimageWeighted(:,:,s);
    RequiredAcc = 1e-10;
    b = sen_map(:,:, s);
    %b=0;
    p = a;
    r = a;
    count = 0;
    delta = r(:)' * r(:) / (a(:)' * a(:));

    while (count < inter_num && delta > RequiredAcc)

        count = count + 1;
            q = (conj(Image) .* Image) .* p;
            if afa ~= 0
                 Lq = Ls * reshape(p, D1 * D2, 1);
                 q = q + afa * reshape(Lq, D1, D2);
            end
          %  Lq=Ls*reshape(q,n1*n2,1);
    %         bs = sum(sum(abs(b)));
    %         if bs ~= 0
    %           % bs=bs/((n1*n2)*1.18^count);
    %            bs = bs / (D1 * D2);
    %         else
    %            bs = 1;
    %         end
            %q=q+afa*reshape(Lq,n1,n2)/bs;

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
    WeightingFunctions(:,:,s) = b;

end
    

