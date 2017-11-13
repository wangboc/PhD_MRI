
function ImgRecon=SEH_GN(trajectory,KspaceData,Density)

%ImgRecon=ifftshift(ifft2(ifftshift(KspaceData.*Density)));
ImgRecon = ifft2(fftshift(KspaceData .* Density));
%ImgRecon=ifft2(fftshift(KspaceData.*Density));