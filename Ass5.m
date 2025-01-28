%cannot be done image deblur c
clear;
% Read the image
img = imread('D:\BENR\YEAR 4 SEM 1\ImgP\assignment\car plate\new\ok list\11.jpg');

% Convert to grayscale if it's a color image
if size(img, 3) == 3
    img = rgb2gray(img);
end

figure (1);
imshow (img);
title ("gray image");


 % Original points 
%P = wnr3; 
% Transformation matrix 
%A = [1.5, 1; 1, 1.5]; % Linear transformation matrix 
%T = [1; 1]; % Translation vector 
% Apply affine transformation 
%P_transformed = (A * P' + T)'; 
%rotated_img = plot([P_transformed(:,1); P_transformed(1,1)], [P_transformed(:,2); P_transformed(1,2)], 'r-o', 'LineWidth', 1.5, 'DisplayName', 'Transformed Shape'); % Rotate by 45 degrees, bilinear 
%figure (4);
%imshow (rotated_img);
%title ("rotated_img");

%figure (4);
%imhist (rotated_img);
%title ("hist");

%I_hist = histeq (wnr3);
%I_hist = medfilt2 (I_hist, [3,3]);


% Sharpened
Isharpen = imsharpen(output_img_y);
figure (2);
imshow (Isharpen);
title ("Sharpened Image");

%crop for the car plate
Icropped = imcrop(Isharpen);
figure (3);
imshow (Icropped);
title ("Cropped Image");

Icropped_resize = imresize (Icropped, 2);
Icropped_filtered = medfilt2 (Icropped_resize, [3,3]);
figure (4);
imshow (Icropped_filtered);
title ("Resized Noise Reduction Image");

%wavelet transformation
% Display the original pixel block 
% Step 1: Row-wise transformation 
[m,n] = size(Icropped_filtered);

if (rem(m,2)) ~= 0
    m = m + 1;
end
if (rem(n,2)) ~= 0
    n = n + 1;
end

Icropped_padded = paddata(Icropped_filtered,[m n]); 
   
LL_R = (Icropped_padded(:, 1:2:end) + Icropped_padded(:, 2:2:end)) / 2; % Approximation 
LH_R = (Icropped_padded(:, 1:2:end) - Icropped_padded(:, 2:2:end)) / 2; % Horizontal Details 

% Step 2: Column-wise transformation 
LL = (LL_R(1:2:end, :) + LL_R(2:2:end, :)) / 2; % Final Approximation (LL) 
LH = (LL_R(1:2:end, :) - LL_R(2:2:end, :)) / 2; % Horizontal Details (LH) 
HL = (LH_R(1:2:end, :) + LH_R(2:2:end, :)) / 2; % Vertical Details (HL) 
HH = (LH_R(1:2:end, :) - LH_R(2:2:end, :)) / 2; % Diagonal Details (HH) 

% Visualization 
figure (5); 
subplot(2, 3, 1); imagesc(Icropped_padded); 
axis equal; 
axis off; colormap gray; 
title('Original Block'); 
subplot(2, 3, 2); imagesc(LL); axis equal; axis off; 
colormap gray; title('LL (Approximation)'); 
subplot(2, 3, 3); imagesc(LH); axis equal; axis off; 
colormap gray; title('LH (Horizontal Details)'); 
subplot(2, 3, 4); imagesc(HL); axis equal; axis off; 
colormap gray; title('HL (Vertical Details)'); 
subplot(2, 3, 5); imagesc(HH); axis equal; axis off; 
colormap gray; title('HH (Diagonal Details)');

%sharpen the car plate
%h = [0 -1 0; -1 5 -1; 0 -1 0];
%Isharpen = imfilter(Icropped_filtered,h);
Isharpen = imsharpen(Icropped_filtered,"Radius",1,"Amount",1,"Threshold",0.23);
%Isharpen = imsharpen(Icropped_filtered);
figure (6);
imshow (Isharpen);
title ("Sharpened Image");

%Image Thresholding
otsu_binary_img = imbinarize(Isharpen, 0.40);
figure (7);
imshow(otsu_binary_img);
title('Image Thresholding');

%Image Inverse
Inversed_image = 1 - otsu_binary_img;
figure (8);
imshow(Inversed_image);
title('Inversed Image');


%ocr detection
ocrResults = ocr(Inversed_image);
ocrResults.Text


%edge detection 
Iedge = edge(otsu_binary_img ,"canny");
figure (9);
imshow (Iedge);
title('Canny Edge Detection Image');

imwrite(LL, 'D:\BENR\YEAR 4 SEM 1\ImgP\assignment\car plate\new\1_compressed.jpg', "Quality", 100) 
imwrite(Icropped_filtered, 'D:\BENR\YEAR 4 SEM 1\ImgP\assignment\car plate\new\1_uncompressed.jpg', "Quality", 100) 
