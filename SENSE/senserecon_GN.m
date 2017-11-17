function recon = senserecon_GN(KspaceDataWeighted, coilprofile, reduc)
%
%   Use coilprofile to form the encoding matrix and find its inverse.
%   Apply the inverse matrix to the aliased image data to recover
%   the original image. Repeat this process for each pixel.
%
%   Input:
%   coilimagealias :  matrices of aliased images
%   coilprofile    :  coil sensitivity maps
%   reduc          :  reduction factor
%
%   Output:
%   recon          :  reconstructed image

coiln = size(coilprofile, 3);
FOV = size(KspaceDataWeighted, 1);

% coilimagealias = zeros(size(KspaceDataWeighted));
for s = 1 : coiln
    coilimagealias(:,:,s) = ifft2(fftshift(KspaceDataWeighted(:, :, s)));
%     观察欠采样后的重叠图像
%     figure; 
%     imshow(rot90(abs(coilimagealias(:, :, s)), -1), []);
% %     观察线圈敏感图像
%     figure; 
%     imshow(rot90(abs(coilprofile(:, :, s)), -1), []);
end

% Img_sos = sqrt(sum(abs(coilimagealias) .^ 2, 3));
% figure; 
% imshow(rot90(abs(Img_sos), -1), []);

recon = zeros(FOV);
for y = 1 : FOV / reduc
    for x = 1 : FOV
        for m = 1 : coiln
            for n = 1 : reduc
                phase = mod((FOV / 2 - FOV / reduc / 2) + y + ...
                    (n - 1) * (FOV / reduc), FOV);
                if phase == 0
                    phase = FOV;
                end
                a(m, 1) = coilimagealias((FOV / 2 - FOV / reduc / 2) + y, x, m);
                encode(m, n) = coilprofile(phase, x, m);
            end
        end

        b = pinv(encode) * a;
        for n = 1 : reduc
            phase = mod((FOV / 2 - FOV / reduc / 2) + y + (n - 1) * (FOV / reduc), FOV);
            if phase == 0
                phase = FOV;
            end
            recon(phase,x) = b(n, 1);
            % imshow(abs(recon), []);
        end
    end
end