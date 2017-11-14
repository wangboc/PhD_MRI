function KspaceTrajectoryData=SEE_GN(trajectory,ImgInput)
%% ¼ÆËã E ÖÐµÄ FT2
% Pruessmann K P, Weiger M, B?rnert P, et al. 
% Advances in sensitivity encoding with arbitrary k©\space trajectories[J]. 
% Magnetic resonance in medicine, 2001, 46(4): 638-651. 
% Figure 1

KspaceTrajectoryData=zeros(size(ImgInput));
kspace_full=fftshift(fft2(ImgInput));
%kspace_full=fftshift(fft2(ImgInput));

KspaceTrajectoryData(trajectory,:)=kspace_full(trajectory,:);



