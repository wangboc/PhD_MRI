function  Fc=PolyfitIFFT_GN(xx, yy, dk1, dk2, R, ACSL, DL, DH, DM, mask, ...
    Order, SI_image, ks_space, Density)

% ************************************************************
%  Matrix
%  这里实际上是对每个线圈获取到的图像进行多项式拟合，并不是敏感度函数
%  因此，当调用estimate_sensitivity_poly_CoilNum函数产生相应WeightingFunctions
%  后，依旧需要计算sum of square
% ************************************************************
[tukey_window,tukey_window_red] = filter_2D(dk2,dk2,R);
[tukey_window2,tukey_window_red2] = filter_img_2D(dk2,dk2,R);

AL = ACSL;%+floor(ACSL/R+0.5);
ks_space_2D = reshape(ks_space, [dk1 + ACSL - DM, dk2]);
%ks_space_2D=ks_space_2D.*tukey_window_red;

for kku = 1 : Order + 1
    for kkv = 1 : Order + 1
        P(1:dk2, 1) = (1:dk2) / dk2 - xx;
        Q(1, 1:dk2) = (1:dk2) / dk2 - yy;
        PQ = P .^ (kku - 1) * Q .^ (kkv - 1);
        kspace = fftshift(fft2(SI_image)) .* tukey_window;
        Img_recon = ifft2(ifftshift(kspace));
        Inimage = PQ .* Img_recon;
                
        kspace_full = fftshift(fft2(Inimage)) .* tukey_window;%.*Density;    
        kspace_red(:, :) = kspace_full(1 : R : dk2, : );
        if AL == 0
           k_space = kspace_red; 
        else
           k_space(1 : floor((DL - 2) / R) + 1, :) = kspace_red(1 : floor((DL - 2) / R) + 1, : );
           k_space(floor((DL - 2) / R) + 2 : floor((DL - 2) / R) + AL + 1, :) = kspace_full(DL : DH, :);
           k_space(floor((DL - 2) / R) + AL + 2 : dk1 + ACSL - DM, :) = ...
               kspace_red(floor((DL - 2) / R) + 2 + DM : dk1, :);
        end        
        dv = (kkv - 1) * (Order + 1) + kku;
        Matrix_E(1 : (dk1 + ACSL - DM) * dk2, dv) = k_space(:);
    end
end
Fc = pinv(Matrix_E) * ks_space_2D(:);
