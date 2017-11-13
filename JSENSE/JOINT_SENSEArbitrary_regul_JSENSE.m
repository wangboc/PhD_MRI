function [ImgRecon,s_poly]=JOINT_SENSEArbitrary_regul_JSENSE(kspace_data,RF,x0,y0,D1,D2,CoilNum,AL,DL,DH,DM,SY,mask,...
                                                         s_poly,order,trajectory,Density,Ls,afa,inter_num,inter)

WeightingFunctions=estimate_sensitivity_poly_CoilNum(D2,D2,CoilNum,x0,y0,s_poly,order); 
[tukey_window,tukey_window_red]=filter_2D(D2,D2,RF);
[tukey_window2,tukey_window_red2]=filter_img_2D(D2,D2,RF);
%figure; mesh(tukey_window2); 

WT_sos=sqrt(sum([abs(WeightingFunctions)].^2,3));
for s = 1 : CoilNum
    Ispace=WeightingFunctions(:,:,s)./WT_sos;
    kspace=fftshift(fft2(Ispace)).*tukey_window;
    WeightingFunctions(:,:,s)=ifft2(ifftshift(kspace));%.*mask+WeightingFunctions0(:,:,s).*(1-mask);    
end
figure;imshow(abs(rot90(WeightingFunctions(:,:,1).*mask,-1)),[]);title('WeightingFunctions of JSENSE');
InitImg=zeros(D2,D2);
ImgRecon=SENSEArbitrary_regul_GN(kspace_data,WeightingFunctions,InitImg,trajectory,Density,Ls,afa,inter_num);
if inter==1
   kspace=fftshift(fft2(ImgRecon)).*tukey_window;%tukey_window2;
else
   kspace=fftshift(fft2(ImgRecon)).*tukey_window;   
end
ImgRecon=ifft2(ifftshift(kspace));

for s = 1 : CoilNum
    Img_WF(:,:,s)=ImgRecon(:,:).*WeightingFunctions(:,:,s);
    kspace_full(:,:,s)=fftshift(fft2(Img_WF(:,:,s))).*tukey_window;%.*Density_WT;
end
SY_VP=ReducedSample(kspace_full,RF,D2,AL,DL,DH,DM);
error_JSENSE=sum(abs(abs(SY)-abs(SY_VP)).^2)/sum(abs(SY).^2)*100
QG=SY-SY_VP;
% *************************************
for s = 1 : CoilNum
    dv = (s-1)*(D1+AL-DM)*D2;
    dk_space(:,1)=SY(dv+1:dv+D2*(D1+AL-DM),1);
    poly_b=PolyfitIFFT_GN(x0,y0,D1,D2,RF,AL,DL,DH,DM,mask,order,ImgRecon,dk_space,Density);  
    s_poly(:,:,s) = reshape(poly_b,[order+1,order+1]);      
end

