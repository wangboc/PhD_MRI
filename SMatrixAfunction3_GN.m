function OutImg = SMatrixAfunction3_GE(InputImg, WeightingFunctions, Sampledtrajectory, I, D)
%calculate the affect of the right I in I*E(H)*D*E*I

InputImg = I .* InputImg;

[D1, D2, nc] = size(WeightingFunctions);
imageweighted = zeros(D1, D2, nc);
xspaceimageWeighted = zeros(D1, D2, nc);

%calculate I*E(H)*D*E
for k = 1 : nc
    imageweighted(:, :, k) = InputImg .* WeightingFunctions(:, :, k);
    kspaceData(:, :, k) = SEE_GN(Sampledtrajectory, imageweighted(:, :, k));
 
    xspaceimage(:,:,k) = SEH_GN(Sampledtrajectory, kspaceData(:, :, k), D);
    xspaceimageWeighted(:, :, k) = xspaceimage(:, :, k) .* conj(WeightingFunctions(:, :, k));
end

OutImg = I .* sum(xspaceimageWeighted, 3); %N^2*1
%OutImg=sum(xspaceimageWeighted,3); %N^2*1

