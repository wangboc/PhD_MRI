
function ImgRecon=SEH_GN(trajectory,KspaceData,Density)
%% ¼ÆËãE(H)ÖÐµÄ FT1¡¢Density correction
% Pruessmann K P, Weiger M, B?rnert P, et al. 
% Advances in sensitivity encoding with arbitrary k©\space trajectories[J]. 
% Magnetic resonance in medicine, 2001, 46(4): 638-651. 
% Figure 1
%%
%ImgRecon=ifftshift(ifft2(ifftshift(KspaceData.*Density)));
ImgRecon = ifft2(fftshift(KspaceData .* Density));
%ImgRecon=ifft2(fftshift(KspaceData.*Density));