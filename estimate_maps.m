function [rough_map, regularization_image]=new_estimate_maps(full_kspace_data,center_lines,...
    center_locs,poly_tick,center_data,Image_size, FEs,PEs);
%Swati Rane, TAMU
%September 1, 2004
%--------------------------------------------------
%---------------------------------------------------
% Copyright (C) 2005 Jim Ji, Texas A&M University.
% All Rights Reserved.

% References:


% Version of 4-June-2005.

% Log:
% Created  2004 Swati Rane	Source code
% Updated
%-----------------------------------------------------
%input= full_kspace_data or reduced kspace_data whatever is available. The
%other is zero
%sampling_locations
%centers of the coil ->center
%variable informing the type of data entered ->data_type
%file containing the extra lines 
%--------------------------------------------------
%output=estimated sensitivity maps->rough_maps
%sum of squares data is nothing but the extra ACs lines used for
%sensitivity estimation-> sos_data
%---------------------------------------------------
%ALL THE PARAMETERS REQUIRED ARE CONTINUOUSLY SAVED IN A FILE CALLED parameters.mat
%THIS FILE IS LOADED BEFORE ANY STEP IS EXECUTED
%---------------------------------------------------
%generate sensitvity maps from the given coil data;
%the center of the k space is collected by filtering with a hamming window
%The unaliased coil image obtained from every coil is then divided by the
%sum os squares image obtained form all the coil images together.
%polynomial filtering is provided,  but the maps look better without that .
%---------------------------------------------------
Ro = FEs;
C = PEs;

k = C / 2;
kc = center_lines / 2; %center_lines/2
w = 20; %20
epsilon = 2;
k1 = center_lines / 2;%center_lines/2  (32-8-15) (24-8-11)  (16-8-7)
option = 2;

%tukey_window = cosine_taper_window(k,kc,w,epsilon,k1,option);
tukey_window = cosine_taper_window(128, 2, 10, 2, 2, 2);
WW = tukey_window(k - center_lines : k + center_lines, k - center_lines : k + center_lines);
ww1 = tukey_window(k : 256, k);


%window_col=hamming(center_lines);
%hamm_window=repmat(window_col,[1 Ro]);
%hamm_window_full=zeros(Ro,C);
%hamm_window_full(:,round(C/2-center_lines/2:C/2+center_lines/2-1))=transpose(hamm_window); %zeroed the figh frequency

% if(center_data~=0)
    [~, ~, num_coils] = size(center_data);
    for i = 1:num_coils
       k_coil(:, :, i) = zeros(Ro, C);
       k_coil(:, center_locs, i) = center_data(:, :, i);
       %k_coil_full(:,:,i)=k_coil(:,:,i).*hamm_window_full;
       k_coil_full(:, :, i) = k_coil(:, :, i) .* tukey_window;
      % figure; mesh(abs(k_coil_full(128-16:128+16,128-16:128+16,i)));
    %   figure; imshow(abs(k_coil_full(128-64:128+64,128-64:128+64,i)),[0,10]);
    %   figure; imshow(abs(k_coil(128-64:128+64,128-64:128+64,i)),[0,10]);
    %   dd
     %  coil_full(:,:,i)=ifftshift(ifft2(ifftshift(k_coil_full(:,:,i))));
       coil_full(:, :, i) = ifft2(ifftshift(k_coil_full(:, :, i)));
    
       %Receiver.sos_data(:,:,i)=center_data(:,:,i);
     end      
%  else
%     [trashp trashf num_coils]=size(full_kspace_data);
%     sampl_in=round(center_lines/Sampling.subsampling_factor);
%     collect_lines=round([C/2-center_lines/2-sampl_in:C/2+center_lines/2+sampl_in-1]);
%     for(i=1:num_coils)
%         k_sub_coil(:,:,i)=zeros(Ro,C);
%         k_sub_coil(:,collect_lines,i)=full_kspace_data(:,collect_lines,i);
%         %Receiver.sos_data(:,:,i)=k_sub_coil(:,collect_lines,i);
%         k_coil_full(:,:,i)=k_sub_coil(:,:,i).*hamm_window_full;
%         coil_full(:,:,i)=ifftshift(ifft2(ifftshift(k_coil_full(:,:,i))));
%     end
%end

%--------sum of squares image--------------------   
%sum_of_squares=sqrt(sum([abs(coil_full)].^2,3)/3);
sum_of_squares = sqrt(sum([abs(coil_full)] .^ 2, 3));
regularization_image = sum_of_squares;
%warning off MATLAB:divideByZero;

for i = 1 : num_coils
    temp = coil_full(:, :, i);
    rough_map(:, :, i) = temp ./ sum_of_squares;
%     figure, imagesc(abs(rough_map(:, :, i))); title('sensitivity estimation with center ACSL')
end

%----------polynomial filtering-------------
if(poly_tick == 1)
    %------------------------
    Np = 1;% order of polynomial
    order = 5; %neighborhood pixels considered for filtering

    %------------------------
    for(i = 1 : num_coils)
        big_map(1 : Ro + 2 * order, 1 : C + 2 * order, i) = 0;
        big_map(order + 1 : order + 1 + Ro - 1, order + 1 : order + 1 + C - 1, i) = rough_map(:, :, i);
        temp_map(:, :) = big_map(:, :, i);
        big_map_contoured(:, :, i) = get_contour(big_map(:, :, i), 0.01) .* big_map(:, :, i);
        clear temp_map;
    end
    roughmap = sensitivity_polyfilt(big_map, Np, order, order);%low pass polynomial
    %for( i=1:num_coils)
    rough_map = roughmap(order + 1 : order + 1 + Ro - 1, order + 1 : order + 1 + C - 1, :);
    %end
end

return