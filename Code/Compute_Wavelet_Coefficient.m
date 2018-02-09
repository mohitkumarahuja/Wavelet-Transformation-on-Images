function Wavelet_Coefficients = Compute_Wavelet_Coefficient(I,J,H0)
% This function is made by Mohit Kumar Ahuja for course Advance Image
% Analysis under the guidance of Dr. Philippe Carre
%
% This function will decompose an Image into J-level wavelet transform
% and will return the Wavelet Coefficients.

%% Computing the High-Pass filter from Low-Pass Filter
for n=0:length(H0)-1  
   G0(n+1)=H0(n+1)*(-1)^(n+1);
end;
G0 = G0(length(G0):-1:1);

for i = 1:J    % The decomposition will run for given J-Level
    
    % The Filter will be applied for Horizontal part as well as for
    % Vertical part so in this code, First we will apply the filter for 
    % horizontal part and then for vertical part.
    
    Image_rows = size(I,1);     % Number of Rows
    Image_columns = size(I,2);  % Number of Columns

    %% Filtering the horizontal part
    
    Horizontal_Part = zeros(size(I)); 
    for i = 1:size(I,1)      % For Rows
        Image_Row = I(i,:);  % One row at a time
        % Low-pass filter
        Image_Row_LPF = pconv(H0,Image_Row);        
        % Downsampling
        Horizontal_Part(i,1:Image_columns/2) = Image_Row_LPF(1:2:length(Image_Row_LPF)); 
        % High-pass filter
        Image_Row_HPF = pconv(G0,Image_Row);
        % Downsampling
        Horizontal_Part(i,Image_columns/2+1:Image_columns)=Image_Row_HPF(1:2:length(Image_Row_HPF));       
    end

    %% Filtering the vertical part
    
    Vertical_Part = zeros(size(I));
    for j = 1:size(I,2)   % For Columns
        Image_Columns = Horizontal_Part(:,j);  % One Column at a time
        % Low-pass Filter
        Image_Columns_LPF = pconv(H0,Image_Columns')'; %Periodic convolution
        % Downsampling
        Vertical_Part(1:Image_rows/2,j) = Image_Columns_LPF(1:2:length(Image_Columns_LPF)); 
        % High-pass Filter
        Image_Columns_HPF = pconv(G0,Image_Columns')';
        % Downsampling
        Vertical_Part(Image_rows/2+1:Image_rows,j) = Image_Columns_HPF(1:2:length(Image_Columns_HPF));
    end
    
    % Making an array of coefficients. Everytime the new coefficients will
    % be concatenated to the previous ones.
    Wavelet_Coefficients(1:Image_rows, 1:Image_columns) = Vertical_Part;
    
    % After every level of decomposition, the Image size reduces to half 
    I = Wavelet_Coefficients(1:Image_rows/2, 1:Image_columns/2);


end
    
end

