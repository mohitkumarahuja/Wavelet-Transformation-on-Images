%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADVANCE IMAGE ANALYSIS
% ----------------------
% Wavelet Transform
% ----------------
% Date: 06 Jan 2018
% Author: Mohit Ahuja !!
% Professor: Dr. Philippe Carre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;


%% PART 1 -> Computing the J-level Wavelet Transform of an NxN Image

% LOAD THE IMAGE
I = double(imread('Cameraman.tif'));  % We have to convert image into double
% I = double(imread('Lena256.bmp'));  % because it will not take Unit8 as an
                                      % input in pconv function.
                                        
% Define WAVELET LEVEL
J = 2;   % Please mention the level of wavelets you want

% Low-pass Filter Analysis
H0 = [0.48296 0.83652 0.22414 -0.12941];

% Function for computing the J-level wavelet transform
Wavelet_Coefficients = Compute_Wavelet_Coefficient(I,J,H0);

% Function for Reconstructing from J-level wavelet transform
Reconstructed_Image = Reconstruct_Using_Wavelet_Coefficient(Wavelet_Coefficients,J,H0);

% Ploting Wavelet Coefficients and Images
figure(1);
subplot(221);
ptr2d(I,0);
colormap(gray(256));
title('Original Image');
subplot(2,2,[2,4]);
ptr2d(Wavelet_Coefficients,2);
title('Wavelet Coefficients');
colormap(gray(256));
subplot(223), affiche(Reconstructed_Image);
title('Reconstructed Image');

% Ploting Error between Original and Reconstructed Image
Error = I - Reconstructed_Image;
figure(2); imshow(Error);
title('Difference between Original and Reconstructed Image');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% END of PART-1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% PART 2 -> Analyse the Denoising Performance of the Wavelet Transform

% Applying Gaussian White Noise
Size = size(I, 1);
Noisy_Image = I+randn(Size,Size)*15;

% Function for computing the J-level wavelet transform
Wavelet_Coefficients_Noise = Compute_Wavelet_Coefficient(Noisy_Image,J,H0);

% Estimate noise level
for S = 1:5
    % Estimation of the noise level, from coeff of first scale
    scale = [Wavelet_Coefficients_Noise(Size/2+1:Size,1:Size/2) Wavelet_Coefficients_Noise(Size/2+1:Size,Size/2+1:Size) Wavelet_Coefficients_Noise(1:Size/2,Size/2+1:Size)];
    sigma = median(abs(scale(:)))/0.85; % Proportion
    threshold = S*sigma;        % level of noise times 3 for the threshold
    % In Gaussian White Noise, 3*std removes all noise
    
    % Apply hard thresholding
    Wavelet_Coefficients_Noise_2 = Wavelet_Coefficients_Noise.*(abs(Wavelet_Coefficients_Noise)>=threshold);
    
    % Apply Soft Thresholding
    % Wavelet_Coefficients_Noise_2 = (sign(Wavelet_Coefficients_Noise).*(abs(Wavelet_Coefficients_Noise)-threshold)).*((abs(Wavelet_Coefficients_Noise)>=threshold));
    
    % Function for Reconstructing from J-level wavelet transform
    Reconstructed_Image_2 = Reconstruct_Using_Wavelet_Coefficient(Wavelet_Coefficients_Noise_2,J,H0);
    
    % Ploting Wavelet Coefficients and Images
    figure(3);
    subplot(221);
    imshow(Noisy_Image,[]);
    title('Original Image');
    subplot(2,2,[2,4]); 
    ptr2d(Wavelet_Coefficients_Noise_2,2); 
    title('Wavelet Coefficients');
    colormap(gray(256));
    subplot(223);
    imshow(Reconstructed_Image_2, []);
    title('Reconstructed Image');
    
    % Compute Peak Signal to Noise Ratio
    [PSNR(S), SNR(S)] = psnr(Reconstructed_Image_2, I);
end

% Ploting SNR with accordance to the number of Standard Deviation
figure(4);
plot(SNR, 'r-*');
title('SNR in accordance to Standard Deviation');
xlabel('No: of STD');
ylabel('SNR (dB)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% END of PART-2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%