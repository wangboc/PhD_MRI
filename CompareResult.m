function [ output ] = CompareResult(Description, Result, sensitivity, StandardImage, kspace_data, ...
    ReduceFactor, D2, ACSL, DL, DH, DM)
%COMPARERESULT �˴���ʾ�йش˺�����ժҪ
%   Description:��Result����
%   Result: �ؽ����
%   StandardImage: ��׼SoS�ο�ͼ��
Img_SoS_NMSE = StandardImage / mean(mean(abs(StandardImage)));
Result_NMSE = Result / mean(mean(abs(Result)));
NMSE = sum(sum(abs(abs(Result_NMSE) - abs(Img_SoS_NMSE)) .^2)) / sum(sum(abs(Img_SoS_NMSE) .^ 2)) * 100;

Img_SoS_SY = ReducedSample(kspace_data, ReduceFactor, D2, ACSL, DL, DH, DM);
Img_WF = zeros(size(Result));
kspace_result = zeros(size(Result));
[~, ~, CoilNum] = size(kspace_data);
for s = 1 : CoilNum
    Img_WF(:, :, s) = Result(:, :) .* sensitivity(:, :, s);
    kspace_result(:, :, s) = fftshift(fft2(Img_WF(:, :, s)));
end
Result_SY = ReducedSample(kspace_result, ReduceFactor, D2, ACSL, DL, DH, DM);
error = sum(abs(abs(Img_SoS_SY) - abs(Result_SY)) .^2) / sum(abs(Img_SoS_SY) .^2) * 100;

output = [Description,':\n', 'error:', num2str(error), '\nNMSE:', num2str(NMSE),'\n**************\n'];
end

