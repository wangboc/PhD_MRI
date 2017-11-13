function  QF0=estimate_kspace_difference(BB20,BB10,WeightingFunctions,sk,dk1,dk2,R,CoilNum,AL,DL,DH,DM,trajectory,area,Ls,afa,inter_num_CG,order)

%  ***********************************************
%  Written by Jinhua Sheng at 09/27/2007 
%  University of Wisconsin
%  **********************************************
QF0=zeros(size(BB10,1)*CoilNum,1);
BB20_Img=ifft2(ifftshift(BB20));
for s = 1 : CoilNum
    if s ~= sk
       B0=(BB20_Img.*WeightingFunctions(:,:,s))./(eps+WeightingFunctions(:,:,sk));
       k_space_full_ekd(:,:,s)=fftshift(fft2(B0));
    end
end
k_space_full_ekd(:,:,sk)=BB20;
KspaceDataWeighted=zeros(dk2,dk2,CoilNum);
KspaceDataWeighted(1:R:dk2,:,:)=k_space_full_ekd(1:R:dk2,:,:);
KspaceDataWeighted(DL:DH,:,:)=k_space_full_ekd(DL:DH,:,:);
InitImg=zeros(dk2,dk2);
Img_space1=JSENSEArbitrary_regul_GN(KspaceDataWeighted,WeightingFunctions,R,InitImg,trajectory,area,Ls,afa,inter_num_CG);

for s = 1 : CoilNum
    Img_WTF1(:,:,s)=Img_space1(:,:).*WeightingFunctions(:,:,s);
    k_space_full_1(:,:,s)=fftshift(fft2(Img_WTF1(:,:,s)));
end
SY1=ReducedSample(k_space_full_1,R,dk2,AL,DL,DH,DM);
QF0((sk-1)*(dk1+AL-DM)*dk2+1:sk*(dk1+AL-DM)*dk2)=BB10(:)-SY1((sk-1)*(dk1+AL-DM)*dk2+1:sk*(dk1+AL-DM)*dk2);
