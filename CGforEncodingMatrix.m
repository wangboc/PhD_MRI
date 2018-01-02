function [ WeightingFunctions ] = CGforEncodingMatrix(kspace_data, Image, mask, Q, pvalue)
% Q �ͷ�ϵ���� pvalue ��ֵ

mask = ~mask;
[m, n, Coilnum] = size(kspace_data);
xspaceimage = zeros(m, n, Coilnum);
xspaceimageWeighted = zeros(m, n, Coilnum);
WeightingFunctions = zeros(m, n, Coilnum);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���²��ù����ݶ��㷨����SENSE���ĸ���C�в鿴
% Pruessmann K P, Weiger M, B?rnert P, et al. 
% Advances in sensitivity encoding with arbitrary k�\space trajectories[J]. 
% Magnetic resonance in medicine, 2001, 46(4): 638-651.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 ������� ��F(H)F)E + F(H)d = 0
% 2 ���ⲿ���гͷ�����ֹ���ֵ��С������ͷ���ֵ
% 2 ����ͷ��� ��F(H)F)E + F(H)d + Q * max{0, -(HE - pvalue)} = 0
% 3 ����� (F(H)F - QH)E + F(H)d + Q * pvalue = 0
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
% %% ������ֵ�����ⲿ��Сֵ����
%     temp = b;
%     b = b .* mask + ~mask;
%     b(abs(b)<pvalue) = pvalue;
%     b = b .* mask + temp .* ~mask;
%% ������
    WeightingFunctions(:,:,s) = b;   
end
    

