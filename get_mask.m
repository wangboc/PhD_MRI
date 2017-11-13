function [mask]=get_mask(I)
% function [mask]=get_contour(I,Threshold)
%  I is the original image (3D or 2D). The salient feature 
%   should be positive.
%  Threshold is for filtering out noise (as 
%  normalized Gaussian). 
%  mask: a binary mask indicate the object/background region
%---------------------------------------------------
% Log:
% Created  2005 Jinhua Sheng	Source code
% Updated
%-----------------------------------------------------

[Km, Km, total_slice]=size(I);

mask = zeros(Km, Km, total_slice);

for index = 1 : total_slice
    Img0 = I(:, :, index);
    Im_max = max(max(abs(Img0(:, :))));
    %[FX,FY] = GRADIENT(Img0);
    [FX, FY] = gradient(Img0);
    FZ = abs(FX) + abs(FY);
    gradient_max = max(max(FZ));
    for kx = 1 : Km
        for ky = 1 : Km
            if FZ(kx, ky) < 0.1 * gradient_max
               Img0_Hist(kx, ky) = 0;
            else
               Img0_Hist(kx, ky) = Img0(kx, ky);
            end 
        end
    end    
    Img0_1D(1 : Km * Km, 1) = Img0_Hist(:);
    nh = hist(Img0_1D, 90);
    %%figure;plot(2:90,abs(nh(2:90)));
    Img_Hist = Img0 >= Im_max * 10 / 90;     %20
    %%figure;imshow(abs(Img_Hist),[]);title('Hist');
    %%figure;imshow(abs(FX)+abs(FY),[]);title('gradient');
    maska = imfill(Img_Hist, 'holes');
    for kx = 1 : Km
        for ky = 1 : Km
            if FZ(kx, ky) >= 0.5 * gradient_max
        %       maska(kx,ky)=1;
            end 
        end
    end    
    maska = bwselect(maska, round(Km / 2), round(Km / 2), 8);
    maska = bwmorph(maska, 'fill',  4);
    maska = bwmorph(maska, 'clean', 4);
    maska = bwmorph(maska, 'close', 7);
end

mask(:, :, index) = maska;
%figure;imshow(FZ/gradient_max-1*mask,[]);title('gradient');
return




