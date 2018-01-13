%% timer on ����ʼ
tic;
%% Initialize parameters and load recon_images
clc; clear all; close all; fclose('all');

path(path,'SENSE'); %Jinhua
path(path,'Sensitivity'); %Jinhua
path(path,'JSENSE'); %Jinhua

global ACSL;
ReduceFactor = 6;  % Facter: R = 1 - 8
ACSL = 64;       %32;  %  for reconstruction
fprintf(['R = ', num2str(ReduceFactor), '       ACLs = ', num2str(ACSL), '\n****************************\n']);
inter_num = 100 ; %2;

global R;
R = ReduceFactor;

inter_num_VP = 6;%6;


order = 10;
FOV = 256;%252;%60;
row = FOV;
column = FOV;
afa = 0.0;%0.0015;

load GE_human_brain;
recon_images = Img;
clear Img;

figure;
imshow(abs(recon_images(:, :, 1)), []);

global D2;
[D1,D2,CoilNum] = size(recon_images);
Img = zeros(D1, D2, CoilNum);
for s = 1 : CoilNum
    Img(:, :, s) = rot90(recon_images(:, :, s), -1);
end
clear recon_images kspace_data
%% �ؽ�SoSͼ������Ӧ����Ȧ���к���
s_poly = zeros(order + 1, order + 1, CoilNum);
Img_SoS = sqrt(sum(abs(Img) .^ 2, 3));

global mask;
[mask] = get_mask(Img_SoS);
Img_NMSE = Img_SoS / mean(mean(abs(Img_SoS)));

figure;
imshow(abs(rot90(Img_SoS, -1)), [0, 0.7 * max(max(abs(Img_SoS)))]);
title('SoS');
WeightingFunctions_standard = zeros(D2, D2, CoilNum);
[tukey_window, tukey_window_red] = filter_2D(D2, D2, ReduceFactor);
% kspace_temp = zeros(D2, D2, CoilNum);

% for s = 1 : CoilNum
%     kspace_temp(:,:,s) = fftshift(fft2(Img(:, :, s)));
% end
% WeightingFunctions_standard = CGforEncodingMatrix(kspace_temp, Img_SoS, 0, 0 ,0, 0);
for s = 1 : CoilNum
    WeightingFunctions_standard(:, :, s) = (Img(:, :, s) ./ (Img_SoS + eps));
end
figure;
imshow(abs(rot90(WeightingFunctions_standard(:, :, 1), 1)), []);
title('Standard WeightingFunctions for Channel 1');

kspace_data = zeros(size(Img));
for s = 1 : CoilNum
    kspace_data(:, :, s) = fftshift(fft2(Img(:, :, s)));
end
%% ����DL��DH��DM�����У�DL��Ӧȫ�������ޣ�DH��Ӧȫ�������ޣ�DM��Ӧȫ����������R����Ӧ��Ƿ����������
Img_reduced = Img(1 : ReduceFactor : D2, :, :);
[D1, D2, CoilNum] = size(Img_reduced);
global DL DH DM;
DL = floor(D2 / 2 + 0.5) + 1 - floor(ACSL / 2 + 0.5);
DH = floor(D2 / 2 + 0.5) + ACSL - floor(ACSL / 2 + 0.5);
if (DL - 1) / R == floor((DL - 1) / R)
    DM = floor((DH - 1) / R) - floor((DL - 1) / R) + 1;
else
    DM = floor((DH - 1) / R) - floor((DL - 1) / R);
end
if ACSL == 0
    DL = floor(D2 / 2 + 0.5) + 1;
    DH = DL;
    DM = 0;
end
%% ����VD Aquisitionͼ�񣬵�Ƶ����ȫ��������Ƶ������� Reduction ȷ���������
trajectory = sort(union(1 : ReduceFactor : D2, DL : DH));
index = trajectory .';
k0 = floor((size(index, 1) + 1) / 2);
area = zeros(size(index, 1), 1);
area(k0) = 1;
for k = 1 : size(index, 1)
    if k > k0
        area(k) = index(k) - index(k - 1);
    end
    if k < k0
        area(k) = index(k + 1) - index(k);
    end
end
Density = zeros(D2,D2);
Density(trajectory, :) = area * ones(1, D2);
figure;
imshow(Density, []);
title(['VD Aquisition with ACSL = ', num2str(ACSL), ' and R = ', num2str(R)]);
% (��δʹ��) JSENSE_tukey_window
%JSENSE_tukey_window = cosine_taper_window(128, 77, 40, 2, 20, 2);
%% estimate maps ͨ������ȫ������Kspace���ݣ�estimate��Ȧ���к�����
center_data = kspace_data(DL : DH, :, :);
[rough_map, regularization_image] = estimate_maps(zeros(D2, D2, CoilNum), ACSL,DL:DH, 0, ...
    permute(center_data, [2 1 3]), 0, D2, D2);
WeightingFunctions = zeros(D2, D2, CoilNum);
for s = 1 : CoilNum
    WeightingFunctions(:, :, s) = permute(rough_map(:, :, s), [2 1]);
end
%% Matrix initialization.ȷ��Kspace����
[~, D2, CoilNum] = size(WeightingFunctions);
kspace_full = kspace_data;


% nencode = ACSL;
% ndim = ceil(D2 / ReduceFactor) * ReduceFactor;
acs_line_location = DL : DH;
KspaceDataWeighted_full = kspace_full;
KspaceDataWeighted = zeros(D2, D2, CoilNum);
KspaceDataWeighted(1 : ReduceFactor : D2, :, :) = KspaceDataWeighted_full(1 : ReduceFactor : D2, :, :);
kspace_full_GN = zeros(size(kspace_full));
for s = 1 : CoilNum
    kspace_full_GN(:, :, s) = kspace_full(:, :, s) .* tukey_window;
end
SY_CG0 = ReducedSample(kspace_full, ReduceFactor, D2, ACSL, DL, DH, DM);
SY = ReducedSample(kspace_full_GN, ReduceFactor, D2, ACSL, DL, DH, DM);
%SY_2D = ReducedSample_2D(kspace_full, ReduceFactor, D2, ACSL, DL, DH, DM);
%% SENSE
InitImg = zeros(D2, D2); 
ImgRecon_SENSE = zeros(D2, D2); 
% ImgRecon_SENSE = senserecon_GN(KspaceDataWeighted, WeightingFunctions, R);
figure,
imshow(abs(rot90(ImgRecon_SENSE, -1)), [0 , 0.7 * max(max(abs(ImgRecon_SENSE)))]);
title('SENSE');
figure;
imshow(rot90(mask, -1), []);
Ls = zeros(0);
if afa ~= 0
    Ls = regularization(D2);
end
figure;
imshow(abs(rot90(WeightingFunctions_standard(:, :, 1) .* mask , -1)), []);
title('SoS�ؽ��ı�׼ WeightingFunctions for Channel 1 with Mask');
figure;
imshow(abs(rot90(WeightingFunctions(:, :, 1) .* mask , -1) ), []);
title('ͨ������ȫ�����������ɵ� WeightingFunctions for Channel 1 with Mask');
%% low resolution reconstruction with CG
kspace_temp = zeros(D2, D2, CoilNum);
kspace_temp(acs_line_location, :, :) = KspaceDataWeighted_full(acs_line_location, :, :);
window = cosine_taper_window(128, 2, 10, 2, 2, 2);
for i = 1:CoilNum
    kspace_temp(:,:,i) = kspace_temp(:,:,i).* window;
end

ImgRecon_CG_low_Resolution = SENSEArbitrary_regul_GN(kspace_temp, WeightingFunctions, InitImg, ...
    trajectory, Density, Ls, afa, inter_num);
ImgRecon_CG_low_Resolution_NMSE = ImgRecon_CG_low_Resolution / mean(mean(abs(ImgRecon_CG_low_Resolution)));
figure,
imshow(abs(rot90(ImgRecon_CG_low_Resolution, -1)), [0, 0.7 * max(max(abs(ImgRecon_CG_low_Resolution)))]);
title('CG (low_resolution)');
%% update mask using ImgRecon_CG_low_Resolution
[mask] = get_mask(ImgRecon_CG_low_Resolution);
% SE=strel('square',3);
% mask=imdilate(mask,SE);
%% Conjugate Gradient for updating sensitivity
Q = 1;        % �ͷ�ϵ��
pvalue = 0.1; % sensitivity ��ֵ����ֹ�����
sensitivity_underCG = CGforEncodingMatrix(kspace_temp, ImgRecon_CG_low_Resolution, 1, mask ,Q, pvalue);
figure,
subplot(2,3,5);
imshow(abs(rot90(WeightingFunctions(:, :, 1), -1)),[]);
title('����ǰ Sensitivity');
subplot(2,3,6);
imshow(abs(rot90(sensitivity_underCG(:, :, 1), -1)),[]);
title(['�ͷ�ϵ��=', num2str(Q), ' ��ֵ=', num2str(pvalue)]);
%% Conjugate Gradient with updated sensitivity
%     for i = 1:CoilNum
%         sensitivity_underCG(:,:,i) = sensitivity_underCG(:,:,i) .* mask;
%     end
KspaceDataWeighted(acs_line_location, :, :) = KspaceDataWeighted_full(acs_line_location, :, :);
ImgRecon_with_updated_sensitivity = SENSEArbitrary_regul_GN(KspaceDataWeighted, sensitivity_underCG, InitImg, ...
trajectory, Density, Ls, afa, inter_num);
% ImgRecon_with_updated_sensitivity = senserecon_GN(KspaceDataWeighted, sensitivity_underCG, R);
subplot(2,3,3);
imshow(abs(rot90(ImgRecon_with_updated_sensitivity.*mask, -1)), []);
title('CG with updated sensitivity');
subplot(2,3,1);
imshow(abs(rot90(Img_SoS.*mask, -1)), []);
title('SoS');
subplot(2,3,4);
imshow(abs(rot90(WeightingFunctions_standard(:,:,1), -1)), []);
title('standard sensitivity');
fprintf(CompareResult('CG with updated sensitivity', ImgRecon_with_updated_sensitivity, sensitivity_underCG, Img_SoS, kspace_data));
%% Standard CG SENSE

ImgRecon_CG = SENSEArbitrary_regul_GN(KspaceDataWeighted, WeightingFunctions, InitImg, ...
trajectory, Density, Ls, afa, inter_num);
subplot(2,3,2);
imshow(abs(rot90(ImgRecon_CG.*mask, -1)), []);
title('CG SENSE');
fprintf(CompareResult('CG SENSE', ImgRecon_CG, WeightingFunctions, Img_SoS, kspace_data));
%% Generate weighting functions and the weighted images. Downsampling by a factor ReduceFactor
% 
return;
k_space_red = zeros(floor(D2 / ReduceFactor), D2, CoilNum);
k_space_red(:, :, :) = kspace_full(1 : ReduceFactor : D2, :, :);
%% WeightingFunctions based on the coefficients of a polynomial
% here xx is the averaged value of x, while it is calculated by
% average(x)/max(x), so does yy.
xx = sum(1 : D2) / D2 ^2; 
yy =    sum(1 : D2) /D2 ^2;
D1 = D2 / R;
for s = 1 : CoilNum
    dv = (s - 1) * (D1 + ACSL - DM) * D2;
    dk_space( :, 1) = SY(dv + 1 : dv + D2 * (D1 + ACSL - DM), 1);
    poly_b = PolyfitIFFT_GN(xx, yy, D1, D2, R, ACSL, DL, DH, DM, ...
        mask, order, ImgRecon_CG, dk_space, Density);
    s_poly(:, :, s) = reshape(poly_b, [order + 1, order + 1]);
end
%% JSENSE
clear kspace_data;
for inter = 1 : inter_num_VP
    inter
    [ImgRecon_JSENSE, s_poly] = JOINT_SENSEArbitrary_regul_JSENSE(KspaceDataWeighted, ...
        R, xx, yy, D1, D2, CoilNum, ACSL, DL, DH, DM, SY, mask, s_poly, order, ...
        trajectory, Density, Ls, afa, inter_num, inter);
    ImgRecon_JSENSE_NMSE = ImgRecon_JSENSE / mean(mean(abs(ImgRecon_JSENSE)));
    JSENSE_NMSE = sum(sum(abs(abs(ImgRecon_JSENSE_NMSE) - abs(Img_NMSE)) .^2)) / ...
        sum(sum(abs(Img_NMSE) .^2)) * 100
    figure,
    imshow(abs(rot90(ImgRecon_JSENSE, -1)), [0, 0.7 * max(max(abs(ImgRecon_JSENSE)))]);
    title(['Iterated Numbers =', num2str(inter), ' using JSENSE']);
end
%% timer off
toc;



