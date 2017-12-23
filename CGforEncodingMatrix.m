function [ WeightingFunctions ] = CGforEncodingMatrix(kspace_data, Image, sen_map, InitImg, trajectory,...
                                        Density, Ls, afa, inter_num, mask, Q)
%CGFORENCODINGMATRIX 此处显示有关此函数的摘要
%   此处显示详细说明

for i = 1 : length(mask(:))
    if(mask(i) == 0)
        mask(i) =  -1;
    end
end
[D1, D2, Coilnum] = size(sen_map);
for s = 1 : Coilnum
    xspaceimage(:, :, s) = ifft2(fftshift(kspace_data(:, :, s)));
    xspaceimageWeighted(:, :, s) =  conj(Image) .* (xspaceimage(:, :, s));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 以下采用共轭梯度算法可在SENSE论文附件C中查看
% Pruessmann K P, Weiger M, B?rnert P, et al. 
% Advances in sensitivity encoding with arbitrary k‐space trajectories[J]. 
% Magnetic resonance in medicine, 2001, 46(4): 638-651.
for s = 1 : Coilnum
    b=sen_map(:,:, s) .* mask;
    a = xspaceimageWeighted(:,:,s) .* mask - Q .* sen_map(:,:, s);
    RequiredAcc = 1e-10;
    p = a;
    r = a;
    count = 0;
    delta = r(:)' * r(:) / (a(:)' * a(:));
    while (count < inter_num && delta > RequiredAcc) 
        count = count + 1;
        q = (conj(Image) .* Image - Q) .* p;
        if afa ~= 0
             Lq = Ls * reshape(p, D1 * D2, 1);
             q = q + afa * reshape(Lq, D1, D2);
        end          
        b = b + (r(:)' * r(:)) / (p(:)' * q(:)) * p;
        tempr = r(:);
        r = r - (r(:)' * r(:)) / (p(:)' * q(:)) * q;
        p = r + (r(:)' * r(:)) / (tempr' * tempr) * p;
        delta = r(:)' * r(:) / (a(:)' * a(:));  
    end
    bb =  b ./ mask;     
    WeightingFunctions(:,:,s) = bb;   
end
    

