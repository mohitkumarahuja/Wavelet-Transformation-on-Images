function Reconstructed_Image = Reconstruct_Using_Wavelet_Coefficient(Wavelet_Coefficients,J,H0)
% This function is made by Mohit Kumar Ahuja for course Advance Image
% Analysis under the guidance of Dr. Philippe Carre
%
% This function will Reconstruct the original Image from J-level wavelet 
% transform and will return the Original Image.

%% Compute high-pass filter
for n=0:length(H0)-1
    G0(n+1)=H0(n+1)*(-1)^(n+1);
end;
G0 = G0(length(G0):-1:1);

% The Reconstrucion filters will be applied to the Horizontal part as well
% as for the Vertical part so in this code, First we will apply the filter
% for vertical part and then for horizontal part.

for N = J-1:-1:0
    Image_rows = size(Wavelet_Coefficients, 1)/(2^N);
    Image_columns = size(Wavelet_Coefficients, 2)/(2^N);
    Wavelet_Coefficients_Level = Wavelet_Coefficients(1:Image_rows, 1:Image_columns);
    
    %% Reconstruction of the Vertical part
    Vertical = zeros(Image_rows, Image_columns);  % Build array of zeros 
    
    for j = 1:Image_columns
        Wave_Coeff_Column = Wavelet_Coefficients_Level(:, j);
        
        % Upsampling
        Wave_Coeff_Column_LPF = zeros(Image_rows, 1); % Build array of zeros
        Wave_Coeff_Column_LPF(1:2:length(Wave_Coeff_Column_LPF)) = Wave_Coeff_Column(1:Image_rows/2); % Fill 1/2 with values from Coefficients
        % Inverse Low-pass Filtering
        Wave_Coeff_Column_LPF = pconv(H0,fliplr(Wave_Coeff_Column_LPF'))'; % fliplr is the reverse process
        % Upsampling
        Wave_Coeff_Column_HPF=zeros(Image_rows, 1);   % Build array of zeros
        Wave_Coeff_Column_HPF(1:2:length(Wave_Coeff_Column_HPF)) = Wave_Coeff_Column(Image_rows/2+1:Image_rows);
        % Inverse High-pass Filtering
        Wave_Coeff_Column_HPF = pconv(G0,fliplr(Wave_Coeff_Column_HPF'))';
        
        % Reconstruct Vertical part
        Vertical(:,j) = fliplr((Wave_Coeff_Column_HPF+Wave_Coeff_Column_LPF)')';
    end
    
    %% Reconstruction of the Horizontal part
    Horizontal = zeros(Image_rows, Image_columns); % Build array of zeros
    
    for i = 1:Image_rows
        Image_Row = Vertical(i,:);
        
        % Upsampling
        Wave_Coeff_Row_LPF = zeros(1, Image_columns);  % Build array of zeros
        Wave_Coeff_Row_LPF(1:2:length(Wave_Coeff_Row_LPF)) = Image_Row(1, 1:Image_columns/2); % Fill 1/2 with values from Coefficients
        % Inverse Lowpass filtering
        Wave_Coeff_Row_LPF = pconv(H0,fliplr(Wave_Coeff_Row_LPF));
        % Upsampling 
        Wave_Coeff_Row_HPF = zeros(1, Image_columns);  % Build array of zeros
        Wave_Coeff_Row_HPF(1:2:length(Wave_Coeff_Row_HPF)) = Image_Row(1, Image_columns/2+1:Image_columns);
        % Inverse Highpass filtering
        Wave_Coeff_Row_HPF = pconv(G0,fliplr(Wave_Coeff_Row_HPF));
        
        % Reconstruct Horizontal part
        Horizontal(i,:) = fliplr((Wave_Coeff_Row_HPF+Wave_Coeff_Row_LPF));
    end
    
    Wavelet_Coefficients(1:Image_rows, 1:Image_columns) = Horizontal;
end

Reconstructed_Image = Wavelet_Coefficients;

end

