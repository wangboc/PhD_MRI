function ImgRec=Sense_GE(k_space_red,WeightingFunctions,ReduceFactor,PiShiftFlag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: SENSE reconstruction using Cartesian samples.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%ImgAliased-The noisy and aliased image in each channel.
%WeightingFunctions-The generated Gaussian-shaped weighting functions.
%ReduceFactor-Reduction factor which determines the sampling interval.
%PiShiftFlag-Determines whether the data set need pi-shifted or not.
%         PiShiftFlag=1 for pi-shift and 0 for no pi-shift.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output:
%ImgRec-The reconstructed image.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 2.0, 02/16/04.
%Written by Dan Xu, Department of Electrical and Computer Engineering,
%University of Illinois at Urbana-Champaign.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get dimensions.
[D1,D2,CoilNum]=size(WeightingFunctions);
for s = 1 : CoilNum
    ImgAliased(:,:,s)= fftshift(ifft(fftshift(k_space_red(:,:,s))));
end
figure, imshow(abs(ImgAliased(:,:,1)),[]);
wqa
%ImgAliased= Pishft1(ifft(Pishft1(k_space_red,1)),1);;
%Solve the superposition equation for each point.
ImgRec=zeros(D1,D2);

for u=1:D2
    for v=1:D1/ReduceFactor
        if PiShiftFlag==0
            WeightingPosition=[v:D1/ReduceFactor:v+D1*(ReduceFactor-1)/ReduceFactor]';
            SenseMatrix=(reshape(WeightingFunctions(WeightingPosition,u,:),length(WeightingPosition),CoilNum)).';
        else
            WeightingPosition=[v+D1/(2*ReduceFactor)+D1/2:D1/ReduceFactor:v+D1*(2*ReduceFactor-1)/(2*ReduceFactor)+D1/2];
            WeightingPosition=1+mod(WeightingPosition-1,D1);
            SenseMatrix=(reshape(WeightingFunctions(WeightingPosition,u,:),length(WeightingPosition),CoilNum)).';
            %SenseMatrix=((-1)^(v+1))*SenseMatrix;
        end
        VectorAliased=ImgAliased(v,u,:);
        VectorAliased=VectorAliased(:);
        %ImgRec(WeightingPosition,u)=pinv(SenseMatrix)*VectorAliased;
        ImgRec(WeightingPosition,u)=((SenseMatrix'*SenseMatrix)\SenseMatrix')*VectorAliased;
     
    end
end