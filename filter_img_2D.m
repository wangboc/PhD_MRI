function  [tukey_window,tukey_window_red]=filter_img_2D(row,colum,R)

% ************************************************************
%  Matrix
% ************************************************************
kc = row / 10; %center_lines/2
w = 20; %100---20
epsilon = 2;
k1 = kc;%center_lines/2  (32-8-15) (24-8-11)  (16-8-7)
option = 2;
tukey_window = cosine_taper_window(row / 2, kc, w, epsilon, k1, option);
%figure; mesh(tukey_window);
Img_space = ifft2(fftshift(tukey_window));
tukey_window = abs(fftshift(fft2(fftshift(Img_space))));
tukey_window_red = tukey_window(1:R:row, :);
