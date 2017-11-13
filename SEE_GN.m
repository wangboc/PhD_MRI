function KspaceTrajectoryData=SEE_GN(trajectory,ImgInput)

KspaceTrajectoryData=zeros(size(ImgInput));
kspace_full=fftshift(fft2(ImgInput));
%kspace_full=fftshift(fft2(ImgInput));

KspaceTrajectoryData(trajectory,:)=kspace_full(trajectory,:);



