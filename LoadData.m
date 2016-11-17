function [I_exact, I_noise] = LoadData(choose)
if choose == 0
    I_exact = imread('cameraman.png');
elseif choose == 1
    I = imread('myphantom1.png');
    I_exact = (1/3)*(I(:,:,1))+(1/3)*(I(:,:,2))+(1/3)*(I(:,:,3));
elseif choose == 2
    I = imread('myphantom2.png');
    I_exact = (1/3)*(I(:,:,1))+(1/3)*(I(:,:,2))+(1/3)*(I(:,:,3));
elseif choose == 3
    I = imread('myphantom3.png');
    I_exact = (1/3)*(I(:,:,1))+(1/3)*(I(:,:,2))+(1/3)*(I(:,:,3));
elseif choose == 4
    I = imread('myphantom4.png');
    I_exact = (1/3)*(I(:,:,1))+(1/3)*(I(:,:,2))+(1/3)*(I(:,:,3));
end
I_exact = double(I_exact);
r = 20;
I_noise = I_exact + r*randn(size(I_exact));
% figure; 
% subplot(1,2,1);
% imshow(I_exact,[]);
% title('Clean image','FontSize',20);
% subplot(1,2,2);
% imshow(I_noise,[]);
% title('Noisy image','FontSize',20);


% % % Construct inpainting region
% % [x,y] = meshgrid(1:size(I_exact,2),1:size(I_exact,1));
% % [th,r] = cart2pol(x-size(I_exact,2)/2,y-size(I_exact,1)/2);
% % D = (sin(r/2+th) > 0.75);   %%%%% the indicator function of the inpainting domain
% % 
% % I_contaminated = I_exact;
% % I_contaminated(D) = rand(nnz(D),1);
% % 
% % subplot(1,3,3);
% % imshow(I_contaminated,[]);
% % title('contaminated','FontSize',20);