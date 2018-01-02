function [ WeightingFunctions ] = CGforEncodingMatrix(kspace_data, Image, mask, Q, pvalue)
% Q 惩罚系数， pvalue 阈值

mask = ~mask;
[m, n, Coilnum] = size(kspace_data);
xspaceimage = zeros(m, n, Coilnum);
xspaceimageWeighted = zeros(m, n, Coilnum);
WeightingFunctions = zeros(m, n, Coilnum);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 以下采用共轭梯度算法可在SENSE论文附件C中查看
% Pruessmann K P, Weiger M, B?rnert P, et al. 
% Advances in sensitivity encoding with arbitrary k‐space trajectories[J]. 
% Magnetic resonance in medicine, 2001, 46(4): 638-651.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 这里求解 （F(H)F)E + F(H)d = 0
% 2 对外部进行惩罚，防止求解值过小，加入惩罚阈值
% 2 加入惩罚项 （F(H)F)E + F(H)d + Q * max{0, -(HE - pvalue)} = 0
% 3 整理后 (F(H)F - QH)E + F(H)d + Q * pvalue = 0
for s = 1 : Coilnum
    xspaceimage(:, :, s) = ifft2(ifftshift(kspace_data(:, :, s)));
    xspaceimageWeighted(:, :, s) =  conj(Image) .* (xspaceimage(:, :, s));
    b = zeros(m,n);
    a = xspaceimageWeighted(:,:,s) + pvalue .* Q .* mask;
    RequiredAcc = 1e-10;
    p = a;
    r = a;
    count = 0;
    delta = r(:)' * r(:) / (a(:)' * a(:));
%     while (count < inter_num && delta > RequiredAcc) 
    while (delta > RequiredAcc) 
        count = count + 1;
        q = (conj(Image) .* Image - Q.* mask) .* p;
        b = b + (r(:)' * r(:)) / (p(:)' * q(:)) * p;
        tempr = r(:);
        r = r - (r(:)' * r(:)) / (p(:)' * q(:)) * q;
        p = r + (r(:)' * r(:)) / (tempr' * tempr) * p;
        delta = r(:)' * r(:) / (a(:)' * a(:));
    end
% %% 根据阈值，对外部过小值修正
%     temp = b;
%     b = b .* mask + ~mask;
%     b(abs(b)<pvalue) = pvalue;
%     b = b .* mask + temp .* ~mask;
%% 结果输出
    WeightingFunctions(:,:,s) = b;   
end
    

