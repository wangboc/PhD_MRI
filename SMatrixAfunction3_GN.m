function OutImg = SMatrixAfunction3_GE(InputImg, WeightingFunctions, Sampledtrajectory, I, D)
%%
% 以下算法可在SENSE论文图1中查看，用于计算 I*E(H)*D*E*I*P
% 可分为三步
% 1：计算 image = I * P
% 2：计算 image` = E(H) * D * E * image
% 3：计算 image`` = I * ∑image`
% Pruessmann K P, Weiger M, B?rnert P, et al. 
% Advances in sensitivity encoding with arbitrary k‐space trajectories[J]. 
% Magnetic resonance in medicine, 2001, 46(4): 638-651.
%% 1、计算 image = I * P
InputImg = I .* InputImg;
%% 2：计算 image` = E(H) * D * E * image
[D1, D2, nc] = size(WeightingFunctions);
imageweighted = zeros(D1, D2, nc);
xspaceimageWeighted = zeros(D1, D2, nc);

for k = 1 : nc
    imageweighted(:, :, k) = InputImg .* WeightingFunctions(:, :, k);
    kspaceData(:, :, k) = SEE_GN(Sampledtrajectory, imageweighted(:, :, k));
 
    xspaceimage(:,:,k) = SEH_GN(Sampledtrajectory, kspaceData(:, :, k), D);
    xspaceimageWeighted(:, :, k) = xspaceimage(:, :, k) .* conj(WeightingFunctions(:, :, k));
end
%% 3：计算 image`` = I * ∑image`
OutImg = I .* sum(xspaceimageWeighted, 3); %N^2*1
