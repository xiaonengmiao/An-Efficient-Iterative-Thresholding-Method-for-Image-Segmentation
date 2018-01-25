close all
clear 

I = imread('../images/cameraman.jpg');
Thres_CV(I,'chan');



I = imread('../images/flowers.jpg'); 
Thres_CV(I,'vector');
 
 

I = imread('../images/flowers.jpg');
Thres_CV(I,'multiphasevector');

 

P = imread('../images/circle.jpg');
% Imnoise the original input
figure()
imshow(P);
I = P;
I(:,:,1) = imnoise(I(:,:,1),'gaussian',0.6,0.5);
% I(:,:,1) = imnoise(I(:,:,1),'salt & pepper',0.8);
Thres_CV(I,'chan');
