callMain=Main();
function [test]=Main()
clear
clc
disp(pwd)
I = imread('D:\project\MAIN PROJECT\Matlab-Project\TestImages\Lena.tiff');
origin_I = double(im2gray(I));
% String to encrypt:
D = double('Very secret string Ashwin J');
D = de2bi(D,8).';
D = D(:).';
num=216
%D=data
%% Set image encryption key and data encryption key
Image_key = 1;
Data_key = 2;
%% Set parameters (convenient for experiment modification)
ref_x = 1; % the number of rows used as reference pixels
ref_y = 1; % the number of columns used as reference pixels
%% Image encryption and data embedding
[encrypt_I,stego_I,emD] = Encrypt_Embed(origin_I,D,Image_key,Data_key,ref_x,ref_y);
%% Data Extraction and Image Recovery
num_emD = length(emD);
if num_emD > 0 % means there is room to embed data
    %-------Extract information in encrypted marked image-------%
    [Side_Information,Refer_Value,Encrypt_exD,Map_I,sign] = Extract_Data(stego_I,num,ref_x,ref_y);
    if sign == 1 % indicates that auxiliary information can be completely extracted
        %---------------Data Decryption----------------%
        [exD] = Encrypt_Data(Encrypt_exD,Data_key);
        % Decrypted string
        rand('seed',Data_key); % set the seed
        E = round(rand(1,numel(Encrypt_exD))*1);
        decrypted = char(sum(reshape(bitxor(Encrypt_exD,E),8,[]).*(2.^(0:7)).'))
        % decrypted = 'Very secret string'
        disp("decrypted")
        disp(decrypted)
%---------------Image Recovery----------------%
        [recover_I] = Recover_Image(stego_I,Image_key,Side_Information,Refer_Value,Map_I,num,ref_x,ref_y);
        %---------------Image comparison----------------%
        figure;
        subplot(221);imshow(origin_I,[]);title('The original image');
        subplot(222);imshow(encrypt_I,[]);title('encrypted image');
        subplot(223);imshow(stego_I,[]);title('Confidential image');
        subplot(224);imshow(recover_I,[]);title('restore image');
        %---------------Result record----------------%
        [m,n] = size(origin_I);
        bpp = num_emD/(m*n);
        %---------------Result judgment----------------%
        check1 = isequal(emD,exD);
        check2 = isequal(origin_I,recover_I);
        if check1 == 1
            disp('Extracting data is exactly the same as embedding data!')
        else
            disp('Warning! Data extraction error!')
        end
        if check2 == 1
            disp('The reconstructed image is exactly the same as the original!')
        else
            disp('Warning! Image reconstruction error!')
        end
        %---------------Result output----------------%
        if check1 == 1 && check2 == 1
            disp(['Embedding capacity equal to : ' num2str(num_emD)])
            disp(['Embedding rate equal to : ' num2str(bpp)])
            fprintf(['The test image------------ OK','\n\n']);
        else
            fprintf(['The test image------------ ERROR','\n\n']);
        end
    else
        disp('Unable to extract all auxiliary information!')
        fprintf(['The test image------------ ERROR','\n\n']);
    end
else
    disp('The auxiliary information is larger than the total embedded amount, so the data cannot be stored!')
    fprintf(['The test image------------ ERROR','\n\n']);
end

test="hi"
end
function value = Binary_Decimalism(bin2_8)
% Convert a binary array to a decimal integer
% Output: bin2_8 (binary array)
% Input: value (decimal integer)
value = 0;
len = length(bin2_8);
for i=1:len
    value = value + bin2_8(i)*(2^(len-i));
end 
end
function [encrypt_I,stego_I,emD] = Encrypt_Embed(origin_I,D,Image_key,Data_key,ref_x,ref_y)
% Function description: Encrypt the original image origin_I and embed the data
% Input: origin_I (original image), D (data to be embedded), Image_key, Data_key (key), ref_x, ref_y (number of rows and columns of reference pixels)
% Output: encrypt_I (encrypted image), stego_I (encrypted marked image), emD (embedded data)
 
% Calculate the predicted value of origin_I
[origin_PV_I] = Predictor_Value(origin_I,ref_x,ref_y);
% label each pixel value (i.e. the position map of the original image)
[Map_origin_I] = Category_Mark(origin_PV_I,origin_I,ref_x,ref_y);
% Huffman-encoded labeling of the labeling categories of pixel values
hist_Map_origin_I = tabulate(Map_origin_I(:)); % Count the number of pixel values ??for each label category
num_Map_origin_I = zeros(9,2);
for i=1:9 % 9 categories of tags
    num_Map_origin_I(i,1) = i-1;
end
[m,~] = size(hist_Map_origin_I);
for i=1:9
    for j=2:m %hist_Map_origin_I The first line counts the number of reference pixels
        if num_Map_origin_I(i,1) == hist_Map_origin_I(j,1)
            num_Map_origin_I(i,2) = hist_Map_origin_I(j,2);
        end
    end
end
[Code,Code_Bin] = Huffman_Code(num_Map_origin_I); % Calculate the mapping relationship of the mark and its binary sequence representation
%% Convert the location map Map_origin_I to a binary array
[Map_Bin] = Map_Binary(Map_origin_I,Code,origin_I,Image_key);
%% Calculate the length of information needed to store the length of Map_Binary
[row,col]=size(origin_I);
max = ceil(log2(row)) + ceil(log2(col)) + 2; % use such a long binary to represent the length of Map_Binary
length_Map = length(Map_Bin);
len_Bin = dec2bin(length_Map)-'0'; % Convert length_Map to binary array
if length(len_Bin) < max
    len = length(len_Bin);
    B = len_Bin;
    len_Bin = zeros(1,max);
    for i=1:len
        len_Bin(max-len+i) = B(i); % Store the length information of Map_Bin
    end
end
%% Auxiliary information required for statistics recovery (Code_Bin, len_Bin, Map_Bin)
Side_Information = [Code_Bin,len_Bin,Map_Bin];
%% encrypt the original image origin_I
[encrypt_I] = Encrypt_Image(origin_I,Image_key);
%% embed information in Encrypt_I
[stego_I,emD] = Embed_Data(encrypt_I,Map_origin_I,Side_Information,D,Data_key,ref_x,ref_y);
end 
function Map_origin_I = Category_Mark(origin_PV_I,origin_I,ref_x,ref_y)
% Function description: mark each pixel value, that is, generate a position map of the original image
% Input: origin_PV_I (predicted value), origin_I (original value), ref_x, ref_y (number of rows and columns of reference pixels)
% Output: Map_origin_I (marker classification of pixel values, i.e. location map)
[row,col] = size(origin_I); % Calculate the row and column values ??of origin_I
Map_origin_I = origin_I; % Build a container that stores the origin_I tag
for i=1:row
    for j=1:col
        if i<=ref_x || j<=ref_y % The first ref_x row and the first ref_y column are used as reference pixels, not marked
            Map_origin_I(i,j) = -1;
        else
            x = origin_I(i,j); % original value
            pv = origin_PV_I(i,j); % predicted value
            for t=7:-1:0
                if floor(x/(2^t)) ~= floor(pv/(2^t))
                    ca = 8-t-1; % marker category used to record pixel values
                    break;
                else
                    ca = 8;
                end
            end
            Map_origin_I(i,j) = ca; % means that the MSB of the ca bit is the same, that is, the ca bit information can be embedded
        end
    end
end
end
function [stego_I,emD] = Embed_Data(encrypt_I,Map_origin_I,Side_Information,D,Data_key,ref_x,ref_y)
% Function description: Embed auxiliary and secret information into encrypted image according to location map
% Input: encrypt_I (encrypted image), Map_origin_I (location map), Side_Information (auxiliary information), D (secret information), Data_key (data encryption key), ref_x, ref_y (number of rows and columns of reference pixels)
% Output: stego_I (encrypted marker image), emD (embedded data)
stego_I = encrypt_I;
[row,col] = size(encrypt_I); % Count the number of rows and columns of encrypt_I
%% Encrypt the original secret information D
[Encrypt_D] = Encrypt_Data(D,Data_key);
%% Record the reference pixels of the first ref_y column and the first ref_x row, and embed them in the image before the secret information
Refer_Value = zeros();
t = 0; % count
for i=1:row
    for j=1:ref_y % before ref_y column
        value = encrypt_I(i,j);
        [bin2_8] = Decimalism_Binary(value); % Convert decimal integer to 8-bit binary array
        Refer_Value(t+1:t+8) = bin2_8;
        t = t + 8;
    end
end
for i=1:ref_x % before ref_x lines
    for j=ref_y+1:col
        value = encrypt_I(i,j);
        [bin2_8] = Decimalism_Binary(value); % Convert decimal integer to 8-bit binary array
        Refer_Value(t+1:t+8) = bin2_8;
        t = t + 8;
    end
end
%% Auxiliary amount
num_D = length(D); % Find the length of the secret information D
num_emD = 0; % count, count the number of embedded secret information
num_S = length(Side_Information); % Find the length of the side information Side_Information
num_side = 0;% count, count the number of embedded auxiliary information
num_RV = length(Refer_Value); % length of reference pixel binary sequence information
num_re = 0; % count, count the length of the binary sequence information embedded in the reference pixel
%% First store auxiliary information in the reference pixels of the first ref_y column and the first ref_x row
for i=1:row
    for j=1:ref_y % before ref_y column
        bin2_8 = Side_Information(num_side+1:num_side+8);
        [value] = Binary_Decimalism(bin2_8); % Convert 8-bit binary array to decimal integer
        stego_I(i,j) = value;
        num_side = num_side + 8;
    end
end
for i=1:ref_x % before ref_x lines
    for j=ref_y+1:col
        bin2_8 = Side_Information(num_side+1:num_side+8);
        [value] = Binary_Decimalism(bin2_8); % Convert 8-bit binary array to decimal integer
        stego_I(i,j) = value;
        num_side = num_side + 8;
    end
end
%% embed auxiliary information, reference pixels and secret data in the remaining positions
for i=ref_x+1:row
    for j=ref_y+1:col
        if num_emD >= num_D % secret data has been embedded
            break;
        end
        %------ means that this pixel can embed 1 bit information------%
        if Map_origin_I(i,j) == 0 %Map=0 means the 1st MSB of the original pixel value is opposite to its predicted value
            if num_side < num_S % auxiliary information is not embedded
                num_side = num_side + 1;
                stego_I(i,j) = mod(stego_I(i,j),2^7) + Side_Information(num_side)*(2^7); % replace 1 MSB
            else
                if num_re < num_RV % reference pixel binary sequence information is not embedded
                    num_re = num_re + 1;
                    stego_I(i,j) = mod(stego_I(i,j),2^7) + Refer_Value(num_re)*(2^7); % replace 1 MSB
                else % embed the secret message at the end
                    num_emD = num_emD + 1;
                    stego_I(i,j) = mod(stego_I(i,j),2^7) + Encrypt_D(num_emD)*(2^7); % replace 1 MSB
                end
            end
%------ means that this pixel can embed 2 bit information------%
        elseif Map_origin_I(i,j) == 1 %Map=1 means the 2 MSB of the original pixel value is opposite to its predicted value
            if num_side < num_S % auxiliary information is not embedded
                if num_side+2 <= num_S % 2 MSBs are used to embed auxiliary information
                    num_side = num_side + 2;
                    stego_I(i,j) = mod(stego_I(i,j),2^6) + Side_Information(num_side-1)*(2^7) + Side_Information(num_side)*(2^6); % replace 2 MSBs
                else
                    num_side = num_side + 1; %1bit auxiliary information
                    num_re = num_re + 1; %1bit reference pixel binary sequence information
                    stego_I(i,j) = mod(stego_I(i,j),2^6) + Side_Information(num_side)*(2^7) + Refer_Value(num_re)*(2^6); % replace 2 MSB
                end
            else
                if num_re < num_RV % reference pixel binary sequence information is not embedded
                    if num_re+2 <= num_RV % 2 MSBs are used to embed reference pixel binary sequence information
                        num_re = num_re + 2;
                        stego_I(i,j) = mod(stego_I(i,j),2^6) + Refer_Value(num_re-1)*(2^7) + Refer_Value(num_re)*(2^6); % replace 2 MSB
                    else
                        num_re = num_re + 1; %1bit reference pixel binary sequence information
                        num_emD = num_emD + 1; %1bit secret information
                        stego_I(i,j) = mod(stego_I(i,j),2^6) + Refer_Value(num_re)*(2^7) + Encrypt_D(num_emD)*(2^6); % replace 2 MSB
                    end
      else
                    if num_emD+2 <= num_D
                        num_emD = num_emD + 2; %2bit secret information
                        stego_I(i,j) = mod(stego_I(i,j),2^6) + Encrypt_D(num_emD-1)*(2^7) + Encrypt_D(num_emD)*(2^6); % replace 2 MSB
                    else
                        num_emD = num_emD + 1; %1bit secret information
                        stego_I(i,j) = mod(stego_I(i,j),2^7) + Encrypt_D(num_emD)*(2^7); % replace 1 MSB
                    end
                end
            end
    %------ means that this pixel can embed 3 bit information------%
        elseif Map_origin_I(i,j) == 2 %Map=2 means that the 3rd MSB of the original pixel value is opposite to its predicted value
            bin2_8 = zeros(1,8); % is used to record the information to be embedded, and the lower bits (LSB) of less than 8 bits default to 0
            if num_side < num_S % auxiliary information is not embedded
                if num_side+3 <= num_S % 3 MSBs are used to embed auxiliary information
                    bin2_8(1:3) = Side_Information(num_side+1:num_side+3);
                    num_side = num_side + 3;
                    [value] = Binary_Decimalism(bin2_8);
                    stego_I(i,j) = mod(stego_I(i,j),2^5) + value; % replace 3 MSBs
                else
                    t = num_S - num_side; % number of remaining auxiliary information
                    bin2_8(1:t) = Side_Information(num_side+1:num_S); %tbit auxiliary information
                    num_side = num_side + t;
                    bin2_8(t+1:3) = Refer_Value(num_re+1:num_re+3-t); %(3-t)bit reference pixel binary sequence information
                    num_re = num_re + 3-t;
                    [value] = Binary_Decimalism(bin2_8);
                    stego_I(i,j) = mod(stego_I(i,j),2^5) + value; % replace 3 MSBs
                end
    else
                if num_re < num_RV % reference pixel binary sequence information is not embedded
                    if num_re+3 <= num_RV % 3 MSBs are used to embed the reference pixel binary sequence information
                        bin2_8(1:3) = Refer_Value(num_re+1:num_re+3);
                        num_re = num_re + 3;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^5) + value; % replace 3 MSBs
                    else
                        t = num_RV - num_re; % number of remaining reference pixel binary sequence information
                        bin2_8(1:t) = Refer_Value(num_re+1:num_RV); %tbit reference pixel binary sequence information
                        num_re = num_re + t;
                        bin2_8(t+1:3) = Encrypt_D(num_emD+1:num_emD+3-t); %(3-t)bit secret information
                        num_emD = num_emD + 3-t;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^5) + value; % replace 3 MSBs
                    end
                        else
                    if num_emD+3 <= num_D
                        bin2_8(1:3) = Encrypt_D(num_emD+1:num_emD+3); %3bit secret information
                        num_emD = num_emD + 3;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^5) + value; % replace 3 MSBs
                    else
                        t = num_D - num_emD; % number of remaining secret messages
                        bin2_8(1:t) = Encrypt_D(num_emD+1:num_emD+t); %tbit secret information
                        num_emD = num_emD + t;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^(8-t)) + value; %Replace t bit MSB
                    end
                end
            end
%------ means that this pixel can embed 4 bit information------%
        elseif Map_origin_I(i,j) == 3 %Map=3 means that the 4th MSB of the original pixel value is opposite to its predicted value
            bin2_8 = zeros(1,8); % is used to record the information to be embedded, and the lower bits (LSB) of less than 8 bits default to 0
            if num_side < num_S % auxiliary information is not embedded
                if num_side+4 <= num_S % 4 MSBs are used to embed auxiliary information
                    bin2_8(1:4) = Side_Information(num_side+1:num_side+4);
                    num_side = num_side + 4;
                    [value] = Binary_Decimalism(bin2_8);
                    stego_I(i,j) = mod(stego_I(i,j),2^4) + value; % replace 4 MSB
                else
                    t = num_S - num_side; % number of remaining auxiliary information
                    bin2_8(1:t) = Side_Information(num_side+1:num_S); %tbit auxiliary information
                    num_side = num_side + t;
                    bin2_8(t+1:4) = Refer_Value(num_re+1:num_re+4-t); %(4-t)bit reference pixel binary sequence information
                    num_re = num_re + 4-t;
                    [value] = Binary_Decimalism(bin2_8);
                    stego_I(i,j) = mod(stego_I(i,j),2^4) + value; % replace 4 MSB
                end
else
                if num_re < num_RV % reference pixel binary sequence information is not embedded
                    if num_re+4 <= num_RV % 4 MSBs are used to embed reference pixel binary sequence information
                        bin2_8(1:4) = Refer_Value(num_re+1:num_re+4);
                        num_re = num_re + 4;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^4) + value; % replace 4 MSB
                    else
                        t = num_RV - num_re; % number of remaining reference pixel binary sequence information
                        bin2_8(1:t) = Refer_Value(num_re+1:num_RV); %tbit reference pixel binary sequence information
                        num_re = num_re + t;
                        bin2_8(t+1:4) = Encrypt_D(num_emD+1:num_emD+4-t); %(4-t)bit secret information
                        num_emD = num_emD + 4-t;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^4) + value; % replace 4 MSB
                    end
        else
                    if num_emD+4 <= num_D
                        bin2_8(1:4) = Encrypt_D(num_emD+1:num_emD+4); %4bit secret information
                        num_emD = num_emD + 4;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^4) + value; % replace 4 MSB
                    else
                        t = num_D - num_emD; % number of remaining secret messages
                        bin2_8(1:t) = Encrypt_D(num_emD+1:num_emD+t); %tbit secret information
                        num_emD = num_emD + t;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^(8-t)) + value; %Replace t bit MSB
                    end
                end
            end 
%------ means that this pixel can embed 5 bit information------%
        elseif Map_origin_I(i,j) == 4 %Map=4 means the 5th MSB of the original pixel value is opposite to its predicted value
            bin2_8 = zeros(1,8); % is used to record the information to be embedded, and the lower bits (LSB) of less than 8 bits default to 0
            if num_side < num_S % auxiliary information is not embedded
                if num_side+5 <= num_S % 5 MSBs are used to embed auxiliary information
                    bin2_8(1:5) = Side_Information(num_side+1:num_side+5);
                    num_side = num_side + 5;
                    [value] = Binary_Decimalism(bin2_8);
                    stego_I(i,j) = mod(stego_I(i,j),2^3) + value; % replace 5 MSB
                else
                    t = num_S - num_side; % number of remaining auxiliary information
                    bin2_8(1:t) = Side_Information(num_side+1:num_S); %tbit auxiliary information
                    num_side = num_side + t;
                    bin2_8(t+1:5) = Refer_Value(num_re+1:num_re+5-t); %(5-t)bit reference pixel binary sequence information
                    num_re = num_re + 5-t;
                    [value] = Binary_Decimalism(bin2_8);
                    stego_I(i,j) = mod(stego_I(i,j),2^3) + value; % replace 5 MSB
                end
             else
                if num_re < num_RV % reference pixel binary sequence information is not embedded
                    if num_re+5 <= num_RV % 5 MSBs are used to embed reference pixel binary sequence information
                        bin2_8(1:5) = Refer_Value(num_re+1:num_re+5);
                        num_re = num_re + 5;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^3) + value; % replace 5 MSB
                    else
                        t = num_RV - num_re; % number of remaining reference pixel binary sequence information
                        bin2_8(1:t) = Refer_Value(num_re+1:num_RV); %tbit reference pixel binary sequence information
                        num_re = num_re + t;
                        bin2_8(t+1:5) = Encrypt_D(num_emD+1:num_emD+5-t); %(5-t)bit secret information
                        num_emD = num_emD + 5-t;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^3) + value; % replace 5 MSB
                    end
                else
                    if num_emD+5 <= num_D
                        bin2_8(1:5) = Encrypt_D(num_emD+1:num_emD+5); %5bit secret information
                        num_emD = num_emD + 5;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^3) + value; % replace 5 MSB
    else
                        t = num_D - num_emD; % number of remaining secret messages
                        bin2_8(1:t) = Encrypt_D(num_emD+1:num_emD+t); %tbit secret information
                        num_emD = num_emD + t;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^(8-t)) + value; %Replace t bit MSB
                    end
                end
            end
    %------ indicates that this pixel can embed 6 bit information------%
        elseif Map_origin_I(i,j) == 5 %Map=5 means the 6MSB of the original pixel value is opposite to its predicted value
            bin2_8 = zeros(1,8); % is used to record the information to be embedded, and the lower bits (LSB) of less than 8 bits default to 0
            if num_side < num_S % auxiliary information is not embedded
                if num_side+6 <= num_S % 6 MSBs are used to embed auxiliary information
                    bin2_8(1:6) = Side_Information(num_side+1:num_side+6);
                    num_side = num_side + 6;
                    [value] = Binary_Decimalism(bin2_8);
                    stego_I(i,j) = mod(stego_I(i,j),2^2) + value; % replace 6-bit MSB
                else
                    t = num_S - num_side; % number of remaining auxiliary information
                    bin2_8(1:t) = Side_Information(num_side+1:num_S); %tbit auxiliary information
                    num_side = num_side + t;
                    bin2_8(t+1:6) = Refer_Value(num_re+1:num_re+6-t); %(6-t)bit reference pixel binary sequence information
                    num_re = num_re + 6-t;
                    [value] = Binary_Decimalism(bin2_8);
                    stego_I(i,j) = mod(stego_I(i,j),2^2) + value; % replace 6-bit MSB
                end
    else
                if num_re < num_RV % reference pixel binary sequence information is not embedded
                    if num_re+6 <= num_RV % 3 MSBs are used to embed the reference pixel binary sequence information
                        bin2_8(1:6) = Refer_Value(num_re+1:num_re+6);
                        num_re = num_re + 6;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^2) + value; % replace 6-bit MSB
                    else
                        t = num_RV - num_re; % number of remaining reference pixel binary sequence information
                        bin2_8(1:t) = Refer_Value(num_re+1:num_RV); %tbit reference pixel binary sequence information
                        num_re = num_re + t;
                        bin2_8(t+1:6) = Encrypt_D(num_emD+1:num_emD+6-t); %(6-t)bit secret information
                        num_emD = num_emD + 6-t;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^2) + value; % replace 6-bit MSB
                    end
                else
                    if num_emD+6 <= num_D
                        bin2_8(1:6) = Encrypt_D(num_emD+1:num_emD+6); %6bit secret information
                        num_emD = num_emD + 6;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^2) + value; % replace 6-bit MSB
    else
                        t = num_D - num_emD; % number of remaining secret messages
                        bin2_8(1:t) = Encrypt_D(num_emD+1:num_emD+t); %tbit secret information
                        num_emD = num_emD + t;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^(8-t)) + value; %Replace t bit MSB
                    end
                end
            end
    %------ means that this pixel can embed 7 bit information------%
        elseif Map_origin_I(i,j) == 6 %Map=6 means the 7th MSB of the original pixel value is opposite to its predicted value
            bin2_8 = zeros(1,8); % is used to record the information to be embedded, and the lower bits (LSB) of less than 8 bits default to 0
            if num_side < num_S % auxiliary information is not embedded
                if num_side+7 <= num_S %7 MSBs are used to embed auxiliary information
                    bin2_8(1:7) = Side_Information(num_side+1:num_side+7);
                    num_side = num_side + 7;
                    [value] = Binary_Decimalism(bin2_8);
                    stego_I(i,j) = mod(stego_I(i,j),2^1) + value; % replace 7 MSB
                else
                    t = num_S - num_side; % number of remaining auxiliary information
                    bin2_8(1:t) = Side_Information(num_side+1:num_S); %tbit auxiliary information
                    num_side = num_side + t;
                    bin2_8(t+1:7) = Refer_Value(num_re+1:num_re+7-t); %(7-t)bit reference pixel binary sequence information
                    num_re = num_re + 7-t;
                    [value] = Binary_Decimalism(bin2_8);
                    stego_I(i,j) = mod(stego_I(i,j),2^1) + value; % replace 7 MSB
                end
    else
                if num_re < num_RV % reference pixel binary sequence information is not embedded
                    if num_re+7 <= num_RV % 7 MSBs are used to embed the reference pixel binary sequence information
                        bin2_8(1:7) = Refer_Value(num_re+1:num_re+7);
                        num_re = num_re + 7;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^1) + value; % replace 7 MSB
                    else
                        t = num_RV - num_re; % number of remaining reference pixel binary sequence information
                        bin2_8(1:t) = Refer_Value(num_re+1:num_RV); %tbit reference pixel binary sequence information
                        num_re = num_re + t;
                        bin2_8(t+1:7) = Encrypt_D(num_emD+1:num_emD+7-t); %(7-t)bit secret information
                        num_emD = num_emD + 7-t;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^1) + value; % replace 7 MSB
                    end
                else
                    if num_emD+7 <= num_D
                        bin2_8(1:7) = Encrypt_D(num_emD+1:num_emD+7); %7bit secret information
                        num_emD = num_emD + 7;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^1) + value; % replace 7 MSB
    else
                        t = num_D - num_emD; % number of remaining secret messages
                        bin2_8(1:t) = Encrypt_D(num_emD+1:num_emD+t); %tbit secret information
                        num_emD = num_emD + t;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^(8-t)) + value; %Replace t bit MSB
                    end
                end
            end
    %------ means that this pixel can embed 8 bit information------%
        elseif Map_origin_I(i,j) == 7 || Map_origin_I(i,j) == 8 %Map=7 means that the 8th MSB (LSB) of the original pixel value is opposite to its predicted value
            bin2_8 = zeros(1,8); % is used to record the information to be embedded, and the lower bits (LSB) of less than 8 bits default to 0
            if num_side < num_S % auxiliary information is not embedded
                if num_side+8 <= num_S %8 MSBs are used to embed auxiliary information
                    bin2_8(1:8) = Side_Information(num_side+1:num_side+8);
                    num_side = num_side + 8;
                    [value] = Binary_Decimalism(bin2_8);
                    stego_I(i,j) = value; % replace 8-bit MSB
                else
                    t = num_S - num_side; % number of remaining auxiliary information
                    bin2_8(1:t) = Side_Information(num_side+1:num_S); %tbit auxiliary information
                    num_side = num_side + t;
                    bin2_8(t+1:8) = Refer_Value(num_re+1:num_re+8-t); %(8-t)bit reference pixel binary sequence information
                    num_re = num_re + 8-t;
                    [value] = Binary_Decimalism(bin2_8);
                    stego_I(i,j) = value; % replace 8-bit MSB
                end
    else
                if num_re < num_RV % reference pixel binary sequence information is not embedded
                    if num_re+8 <= num_RV % 8 MSBs are used to embed reference pixel binary sequence information
                        bin2_8(1:8) = Refer_Value(num_re+1:num_re+8);
                        num_re = num_re + 8;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = value; % replace 8-bit MSB
                    else
                        t = num_RV - num_re; % number of remaining reference pixel binary sequence information
                        bin2_8(1:t) = Refer_Value(num_re+1:num_RV); %tbit reference pixel binary sequence information
                        num_re = num_re + t;
                        bin2_8(t+1:8) = Encrypt_D(num_emD+1:num_emD+8-t); %(8-t)bit secret information
                        num_emD = num_emD + 8-t;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = value; % replace 8-bit MSB
                    end
    else
                    if num_emD+8 <= num_D
                        bin2_8(1:8) = Encrypt_D(num_emD+1:num_emD+8); %8bit secret information
                        num_emD = num_emD + 8;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = value; % replace 8-bit MSB
                    else
                        t = num_D - num_emD; % number of remaining secret messages
                        bin2_8(1:t) = Encrypt_D(num_emD+1:num_emD+t); %tbit secret information
                        num_emD = num_emD + t;
                        [value] = Binary_Decimalism(bin2_8);
                        stego_I(i,j) = mod(stego_I(i,j),2^(8-t)) + value; %Replace t bit MSB
                    end
                end
            end
        end
    end
end
%% Statistics embedded secret data
emD = D(1:num_emD);
end
function [Encrypt_D] = Encrypt_Data(D,Data_key)
% Encryption binary array
rand('seed',Data_key); % set the seed
E = round(rand(1,numel(D))*1);
Encrypt_D=D
% crypted string
 Encrypt_D = bitxor(D,E);
disp(length(Encrypt_D))
end
function [Code,Code_Bin] = Huffman_Code(num_Map_origin_I)
% Function description: use variable-length encoding (multi-bit 0/1 encoding) to represent the marker category of pixel values
% Input: num_Map_origin_I (statistics of pixel value marker categories)
% Output: Code (mapping relationship), Code_Bin (binary representation of Code)
% Remarks: Use the 9 codes of {00, 01, 100, 101, 1100, 1101, 1110, 11110, 11111} to represent 9 types of markers
% Rule: The more pixels in the marked category, the shorter the encoding length used to represent its category
%{00,01,100,101,1100,1101,1110,11110,11111}?{0,1,4,5,12,13,14,30,31}
%% Find its mapping encoding relationship
Code = [-1,0;-1,1;-1,4;-1,5;-1,12;-1,13;-1,14;-1,30;-1,31];
for i=1:9
    drder=1;
    for j=1:9
        if num_Map_origin_I(i,2) < num_Map_origin_I(j,2)
            drder = drder + 1;
        end
    end
    while Code(drder) ~= -1 % prevent the number of pixels in the two marker classes from being equal
        drder = drder + 1;
    end
    Code(drder,1) = num_Map_origin_I(i,1);
end
%% Represent the Map mapping relationship as a binary bit stream
Code_Bin = zeros();
t = 0; % count
for i=0:8
    for j=1:9
        if Code(j,1) == i
            value = Code(j,2);
        end
    end
    if value == 0
        Code_Bin(t+1) = 0;
        Code_Bin(t+2) = 0;
        t = t+2;
    elseif value == 1
        Code_Bin(t+1) = 0;
        Code_Bin(t+2) = 1;
        t = t+2;
    else
        add = ceil(log2(value+1)); % indicates the length of the marker encoding
        Code_Bin(t+1:t+add) = dec2bin(value)-'0'; % Convert value to binary array
        t = t + add;
    end
end
end
function [Map_Bin] = Map_Binary(Map_origin_I,Code,origin_I,Image_key)
% Function description: Convert the location map Map_origin_I into a binary array Map
% Input: Map_origin_I (position map of the original image), Code (mapping relationship)
% Output: Map_Bin (binary array of original image location map)
[row,col] = size(Map_origin_I); % Calculate the row and column values ??of Map_origin_II
Map_Bin = zeros();
t = 0; % count, length of binary array
for i=1:row
    for j=1:col
        if Map_origin_I(i,j) == -1 % The point marked with -1 is used as a reference pixel and is not counted
            continue;
        end
        for k=1:9
            if Map_origin_I(i,j) == Code(k,1)
                value = Code(k,2);
                break;
            end
        end
        if value == 0
            Map_Bin(t+1) = 0;
            Map_Bin(t+2) = 0;
            t = t+2;
        elseif value == 1
            Map_Bin(t+1) = 0;
            Map_Bin(t+2) = 1;
            t = t+2;
        else
            add = ceil(log2(value+1)); % indicates the length of the marker encoding
            Map_Bin(t+1:t+add) = dec2bin(value)-'0'; % Convert value to binary array
            t = t + add;
        end
    end
end
 
% Function description: Bit-level XOR encryption of image origin_I
% Input: origin_I (original image), Image_key (image encryption key)
% output: encrypt_I (encrypted image)
[row,col] = size(origin_I); % Calculate the row and column values ??of origin_I
encrypt_I = origin_I; % build a container that stores encrypted images
%% Generate a random matrix of the same size as origin_I based on the key
rand('seed',Image_key); %Set the seed
E = round(rand(row,col)*255); % Randomly generate row*col matrix
%% Bit-level encryption of image origin_I according to E
for i=1:row
    for j=1:col
        encrypt_I(i,j) = bitxor(origin_I(i,j),E(i,j));
    end
end
end
function [encrypt_I] = Encrypt_Image(origin_I,Image_key)
% Function description: Bit-level XOR encryption of image origin_I
% Input: origin_I (original image), Image_key (image encryption key)
% output: encrypt_I (encrypted image)
[row,col] = size(origin_I); % Calculate the row and column values ??of origin_I
encrypt_I = origin_I; % build a container that stores encrypted images
%% Generate a random matrix of the same size as origin_I based on the key
rand('seed',Image_key); %Set the seed
E = round(rand(row,col)*255); % Randomly generate row*col matrix
%% Bit-level encryption of image origin_I according to E
for i=1:row
    for j=1:col
        encrypt_I(i,j) = bitxor(origin_I(i,j),E(i,j));
    end
end
end
function origin_PV_I = Predictor_Value(origin_I,ref_x,ref_y)
           disp("In the function") 
           [row,col] = size(origin_I); 
            origin_PV_I = origin_I;  
            for i=ref_x+1:row  
                for j=ref_y+1:col  
                     a = origin_I(i-1,j);
                     b = origin_I(i-1,j-1);
                     c = origin_I(i,j-1);
                    if b <= min(a,c)
                        origin_PV_I(i,j) = max(a,c);
                    elseif b >= max(a,c)
                        origin_PV_I(i,j) = min(a,c);
                    else
                        origin_PV_I(i,j) = a + c - b;
                    end
               end 
            end
            disp("out of the function") 
            
end
function [bin2_8] = Decimalism_Binary(value)
% º¯ÊýËµÃ÷£º½«Ê®½øÖÆ»Ò¶ÈÏñËØÖµ×ª»»³É8Î»¶þ½øÖÆÊý×é
% ÊäÈë£ºvalue£¨Ê®½øÖÆ»Ò¶ÈÏñËØÖµ£©
% Êä³ö£ºbin2_8£¨8Î»¶þ½øÖÆÊý×é£©
bin2_8 = dec2bin(value)-'0';
if length(bin2_8) < 8
    len = length(bin2_8);
    B = bin2_8;
    bin2_8 = zeros(1,8);
    for i=1:len
        bin2_8(8-len+i) = B(i); %²»×ã8Î»Ç°Ãæ²¹³ä0
    end 
end
 end
function [Side_Information,Refer_Value,Encrypt_exD,Map_I,sign] = Extract_Data(stego_I,num,ref_x,ref_y)
% Function Description: Extract information from encrypted labeled images
% Input: stego_I (encrypted marked image), num (length of secret information), ref_x, ref_y (number of rows and columns of reference pixels)
% Output: Side_Information (auxiliary information), Refer_Value (reference pixel information), Encrypt_exD (encrypted secret information), Map_I (location map), sign (judgment mark)
[row,col]=size(stego_I); % count the number of rows and columns of stego_I
%% Build a matrix to store the location map
Map_I = zeros(row,col); % Build a matrix that stores the location map
for i=1:row
    for j=1:ref_y
        Map_I(i,j) = -1; % The previous ref_y column is the reference pixel, and no marking is performed
    end
end
for i=1:ref_x
    for j=ref_y+1:col
        Map_I(i,j) = -1; % The previous ref_x line is the reference pixel, not marked
    end
end
%% First extract the auxiliary information in the first ref_y column and the first ref_x row
Side_Information = zeros();
num_side = 0;% count, count the number of auxiliary information extracted
for i=1:row
    for j=1:ref_y
        value = stego_I(i,j);
        [bin2_8] = Decimalism_Binary(value); % Convert decimal integer to 8-bit binary array
        Side_Information(num_side+1:num_side+8) = bin2_8;
        num_side = num_side + 8;
    end
end
for i=1:ref_x
    for j=ref_y+1:col
        value = stego_I(i,j);
        [bin2_8] = Decimalism_Binary(value); % Convert decimal integer to 8-bit binary array
        Side_Information(num_side+1:num_side+8) = bin2_8;
        num_side = num_side + 8;
    end
end
%% Extract auxiliary information representing mapping rules
Code_Bin = Side_Information(1:32); % The first 32 bits are mapping rule information
Code = [0,-1;1,-1;2,-1;3,-1;4,-1;5,-1;6,-1;7,-1;8,-1];
this_end = 0;
for i=1:9 % Convert binary sequence map to integer map
    last_end = this_end;
    [code_value,this_end] = Huffman_DeCode(Code_Bin,last_end);
    Code(i,2) = code_value;
end
%% Extract the length information of the binary sequence of the position map
max = ceil(log2(row)) + ceil(log2(col)) + 2; % use such a long binary representation Map_I to convert the length of the binary sequence
len_Bin = Side_Information(33:32+max); % The first 33 to 32+max bits are the length information of the binary sequence of the position map
num_Map = 0; % Convert the binary sequence len_Bin to a decimal number
for i=1:max
    num_Map = num_Map + len_Bin(i)*(2^(max-i));
end
%% Auxiliary amount
num_S = 32 + max + num_Map; % length of auxiliary information
Refer_Value = zeros();
num_RV = (ref_x*row+ref_y*col-ref_x*ref_y)*8; % length of reference pixel binary sequence information
num_re = 0; % count, count and extract the length of the binary sequence information of the reference pixel
Encrypt_exD = zeros();
num_D = num; % length of binary secret message
num_exD = 0; % count, count the number of embedded secret information
%% Extract information beyond the first rows and columns
this_end = 32 + max; % the auxiliary information in front of is not a location map
sign = 1; % indicates that the data can be completely extracted to restore the image
for i=ref_x+1:row
    if sign == 0 % means that the data cannot be fully extracted to restore the image
        break;
    end
    for j=ref_y+1:col
        if num_exD >= num_D % secret data has been extracted
            break;
        end
        %------ Convert the current decimal pixel value to an 8-bit binary array------%
        value = stego_I(i,j);
        [bin2_8] = Decimalism_Binary(value);
        %--Calculate how many bits of information the current pixel can extract through auxiliary information--%
        last_end = this_end;
        [map_value,this_end] = Huffman_DeCode(Side_Information,last_end);
        if map_value == -1 % indicates that the length of auxiliary information is not enough to restore the next Huffman code
sign = 0;
            break;
        end
        for k=1:9
            if map_value == Code(k,2)
                Map_I(i,j) = Code(k,1); % Position map information of the current pixel
                break;
            end
        end
        %------- means that this pixel can extract 1 bit information-------%
        if Map_I(i,j) == 0 %Map=0 means that the 1st MSB of the original pixel value is opposite to its predicted value
            if num_side < num_S % auxiliary information has not been extracted
                num_side = num_side + 1;
                Side_Information(num_side) = bin2_8(1);
            else
                if num_re < num_RV % reference pixel binary sequence information has not been extracted
                    num_re = num_re + 1;
                    Refer_Value(num_re) = bin2_8(1);
                else % Finally extract the secret information
                    num_exD = num_exD + 1;
                    Encrypt_exD(num_exD) = bin2_8(1);
                end
end
        %------- means that this pixel can extract 2 bit information-------%
        elseif Map_I(i,j) == 1 %Map=1 means the 2 MSB of the original pixel value is opposite to its predicted value
            if num_side < num_S % auxiliary information has not been extracted
                if num_side+2 <= num_S %2 MSBs are auxiliary information
                    Side_Information(num_side+1:num_side+2) = bin2_8(1:2);
                    num_side = num_side + 2;
                else
                    num_side = num_side + 1; %1bit auxiliary information
                    Side_Information(num_side) = bin2_8(1);
                    num_re = num_re + 1; %1bit reference pixel binary sequence information
                    Refer_Value(num_re) = bin2_8(2);
                end
            else
                if num_re < num_RV % reference pixel binary sequence information has not been extracted
                    if num_re+2 <= num_RV % 2 MSBs are reference pixel binary sequence information
                        Refer_Value(num_re+1:num_re+2) = bin2_8(1:2);
                        num_re = num_re + 2;
                    else
                        num_re = num_re + 1; %1bit reference pixel binary sequence information
Refer_Value(num_re) = bin2_8(1);
                        num_exD = num_exD + 1; %1bit secret information
                        Encrypt_exD(num_exD) = bin2_8(2);
                    end
                else
                    if num_exD+2 <= num_D
                        Encrypt_exD(num_exD+1:num_exD+2) = bin2_8(1:2); %2bit secret information
                        num_exD = num_exD + 2;
                    else
                        num_exD = num_exD + 1; %1bit secret information
                        Encrypt_exD(num_exD) = bin2_8(1);
                    end
                end
            end
        %------- means that this pixel can extract 3 bit information-------%
        elseif Map_I(i,j) == 2 %Map=2 means that the 3rd MSB of the original pixel value is opposite to its predicted value
            if num_side < num_S % auxiliary information has not been extracted
                if num_side+3 <= num_S % 3 MSBs are auxiliary information
                    Side_Information(num_side+1:num_side+3) = bin2_8(1:3);
                    num_side = num_side + 3;
                else
                    t = num_S - num_side; % number of remaining auxiliary information
                    Side_Information(num_side+1:num_side+t) = bin2_8(1:t); %tbit auxiliary information
                    num_side = num_side + t;
                    Refer_Value(num_re+1:num_re+3-t) = bin2_8(t+1:3); %(3-t)bit reference pixel binary sequence information
                    num_re = num_re + 3-t;
                end
            else
                if num_re < num_RV % reference pixel binary sequence information has not been extracted
if num_re+3 <= num_RV % 3 MSBs are reference pixel binary sequence information
                        Refer_Value(num_re+1:num_re+3) = bin2_8(1:3);
                        num_re = num_re + 3;
                    else
                        t = num_RV - num_re; % number of remaining reference pixel binary sequence information
                        Refer_Value(num_re+1:num_re+t) = bin2_8(1:t); %tbit reference pixel binary sequence information
                        num_re = num_re + t;
                        Encrypt_exD(num_exD+1:num_exD+3-t) = bin2_8(t+1:3); %(3-t)bit secret information
                        num_exD = num_exD + 3-t;
                    end
                else
                    if num_exD+3 <= num_D
                        Encrypt_exD(num_exD+1:num_exD+3) = bin2_8(1:3); %3bit secret information
                        num_exD = num_exD + 3;
                    else
                        t = num_D - num_exD;
                        Encrypt_exD(num_exD+1:num_exD+t) = bin2_8(1:t); %tbit secret information
                        num_exD = num_exD + t;
                    end
                end
            end
        %------- means that this pixel can extract 4 bit information-------%
        elseif Map_I(i,j) == 3 %Map=3 means the 4th MSB of the original pixel value is opposite to its predicted value
            if num_side < num_S % auxiliary information has not been extracted
                if num_side+4 <= num_S %4 MSBs are auxiliary information
                    Side_Information(num_side+1:num_side+4) = bin2_8(1:4);
                    num_side = num_side + 4;
                else
t = num_S - num_side; % number of remaining auxiliary information
                    Side_Information(num_side+1:num_side+t) = bin2_8(1:t); %tbit auxiliary information
                    num_side = num_side + t;
                    Refer_Value(num_re+1:num_re+4-t) = bin2_8(t+1:4); %(4-t)bit reference pixel binary sequence information
                    num_re = num_re + 4-t;
                end
            else
                if num_re < num_RV % reference pixel binary sequence information has not been extracted
                    if num_re+4 <= num_RV % 4 MSBs are reference pixel binary sequence information
                        Refer_Value(num_re+1:num_re+4) = bin2_8(1:4);
                        num_re = num_re + 4;
                    else
                        t = num_RV - num_re; % number of remaining reference pixel binary sequence information
                        Refer_Value(num_re+1:num_re+t) = bin2_8(1:t); %tbit reference pixel binary sequence information
                        num_re = num_re + t;
                        Encrypt_exD(num_exD+1:num_exD+4-t) = bin2_8(t+1:4); %(4-t)bit secret information
                        num_exD = num_exD + 4-t;
                    end
                else
                    if num_exD+4 <= num_D
                        Encrypt_exD(num_exD+1:num_exD+4) = bin2_8(1:4); %4bit secret information
                        num_exD = num_exD + 4;
                    else
                        t = num_D - num_exD;
                        Encrypt_exD(num_exD+1:num_exD+t) = bin2_8(1:t); %tbit secret information
                        num_exD = num_exD + t;
                    end
                end
            end
%------- means that this pixel can extract 5 bit information-------%
        elseif Map_I(i,j) == 4 %Map=4 means that the 5th MSB of the original pixel value is the opposite of its predicted value
            if num_side < num_S % auxiliary information has not been extracted
                if num_side+5 <= num_S %5 MSBs are auxiliary information
                    Side_Information(num_side+1:num_side+5) = bin2_8(1:5);
                    num_side = num_side + 5;
                else
                    t = num_S - num_side; % number of remaining auxiliary information
                    Side_Information(num_side+1:num_side+t) = bin2_8(1:t); %tbit auxiliary information
                    num_side = num_side + t;
                    Refer_Value(num_re+1:num_re+5-t) = bin2_8(t+1:5); %(5-t)bit reference pixel binary sequence information
                    num_re = num_re + 5-t;
                end
            else
                if num_re < num_RV % reference pixel binary sequence information has not been extracted
                    if num_re+5 <= num_RV % 5 bits MSB are reference pixel binary sequence information
                        Refer_Value(num_re+1:num_re+5) = bin2_8(1:5);
                        num_re = num_re + 5;
                    else
                        t = num_RV - num_re; % number of remaining reference pixel binary sequence information
                        Refer_Value(num_re+1:num_re+t) = bin2_8(1:t); %tbit reference pixel binary sequence information
                        num_re = num_re + t;
                        Encrypt_exD(num_exD+1:num_exD+5-t) = bin2_8(t+1:5); %(5-t)bit secret information
                        num_exD = num_exD + 5-t;
                    end
else
                    if num_exD+5 <= num_D
                        Encrypt_exD(num_exD+1:num_exD+5) = bin2_8(1:5); %5bit secret information
                        num_exD = num_exD + 5;
                    else
                        t = num_D - num_exD;
                        Encrypt_exD(num_exD+1:num_exD+t) = bin2_8(1:t); %tbit secret information
                        num_exD = num_exD + t;
                    end
                end
            end
            %------- means that this pixel can extract 6 bit information-------%
        elseif Map_I(i,j) == 5 %Map=5 means the 6MSB of the original pixel value is opposite to its predicted value
            if num_side < num_S % auxiliary information has not been extracted
                if num_side+6 <= num_S % 6 MSBs are auxiliary information
                    Side_Information(num_side+1:num_side+6) = bin2_8(1:6);
                    num_side = num_side + 6;
                else
                    t = num_S - num_side; % number of remaining auxiliary information
                    Side_Information(num_side+1:num_side+t) = bin2_8(1:t); %tbit auxiliary information
                    num_side = num_side + t;
                    Refer_Value(num_re+1:num_re+6-t) = bin2_8(t+1:6); %(6-t)bit reference pixel binary sequence information
                    num_re = num_re + 6-t;
                end
else
                if num_re < num_RV % reference pixel binary sequence information has not been extracted
                    if num_re+6 <= num_RV % 6 MSBs are reference pixel binary sequence information
                        Refer_Value(num_re+1:num_re+6) = bin2_8(1:6);
                        num_re = num_re + 6;
                    else
                        t = num_RV - num_re; % number of remaining reference pixel binary sequence information
                        Refer_Value(num_re+1:num_re+t) = bin2_8(1:t); %tbit reference pixel binary sequence information
                        num_re = num_re + t;
                        Encrypt_exD(num_exD+1:num_exD+6-t) = bin2_8(t+1:6); %(6-t)bit secret information
                        num_exD = num_exD + 6-t;
                    end
                else
                    if num_exD+6 <= num_D
                        Encrypt_exD(num_exD+1:num_exD+6) = bin2_8(1:6); %6bit secret information
                        num_exD = num_exD + 6;
                    else
                        t = num_D - num_exD;
                        Encrypt_exD(num_exD+1:num_exD+t) = bin2_8(1:t); %tbit secret information
                        num_exD = num_exD + t;
                    end
end
            end
            %------- means that this pixel can extract 7 bit information-------%
        elseif Map_I(i,j) == 6 %Map=6 means the 7th MSB of the original pixel value is the opposite of its predicted value
            if num_side < num_S % auxiliary information has not been extracted
                if num_side+7 <= num_S %7 MSBs are auxiliary information
                    Side_Information(num_side+1:num_side+7) = bin2_8(1:7);
                    num_side = num_side + 7;
                else
                    t = num_S - num_side; % number of remaining auxiliary information
                    Side_Information(num_side+1:num_side+t) = bin2_8(1:t); %tbit auxiliary information
                    num_side = num_side + t;
                    Refer_Value(num_re+1:num_re+7-t) = bin2_8(t+1:7); %(7-t)bit reference pixel binary sequence information
                    num_re = num_re + 7-t;
                end
            else
                if num_re < num_RV % reference pixel binary sequence information has not been extracted
                    if num_re+7 <= num_RV % 7 bits MSB are reference pixel binary sequence information
                        Refer_Value(num_re+1:num_re+7) = bin2_8(1:7);
                        num_re = num_re + 7;
                    else
                        t = num_RV - num_re; % number of remaining reference pixel binary sequence information
                        Refer_Value(num_re+1:num_re+t) = bin2_8(1:t); %tbit reference pixel binary sequence information
                        num_re = num_re + t;
Encrypt_exD(num_exD+1:num_exD+7-t) = bin2_8(t+1:7); %(7-t)bit secret information
                        num_exD = num_exD + 7-t;
                    end
                else
                    if num_exD+7 <= num_D
                        Encrypt_exD(num_exD+1:num_exD+7) = bin2_8(1:7); %7bit secret information
                        num_exD = num_exD + 7;
                    else
                        t = num_D - num_exD;
                        Encrypt_exD(num_exD+1:num_exD+t) = bin2_8(1:t); %tbit secret information
                        num_exD = num_exD + t;
                    end
                end
            end
            %------- means that this pixel can extract 8 bit information-------%
        elseif Map_I(i,j) == 7 || Map_I(i,j) == 8 %Map=7 means that the 8th MSB (LSB) of the original pixel value is opposite to its predicted value
            if num_side < num_S % auxiliary information has not been extracted
                if num_side+8 <= num_S %8 MSBs are auxiliary information
                    Side_Information(num_side+1:num_side+8) = bin2_8(1:8);
                    num_side = num_side + 8;
                else
                    t = num_S - num_side; % number of remaining auxiliary information
                    Side_Information(num_side+1:num_side+t) = bin2_8(1:t); %tbit auxiliary information
                    num_side = num_side + t;
                    Refer_Value(num_re+1:num_re+8-t) = bin2_8(t+1:8); %(8-t)bit reference pixel binary sequence information
                    num_re = num_re + 8-t;
                end
else
                if num_re < num_RV % reference pixel binary sequence information has not been extracted
                    if num_re+8 <= num_RV % 8 MSBs are reference pixel binary sequence information
                        Refer_Value(num_re+1:num_re+8) = bin2_8(1:8);
                        num_re = num_re + 8;
                    else
                        t = num_RV - num_re; % number of remaining reference pixel binary sequence information
                        Refer_Value(num_re+1:num_re+t) = bin2_8(1:t); %tbit reference pixel binary sequence information
                        num_re = num_re + t;
                        Encrypt_exD(num_exD+1:num_exD+8-t) = bin2_8(t+1:8); %(8-t)bit secret information
                        num_exD = num_exD + 8-t;
                    end
                else
                    if num_exD+8 <= num_D
                        Encrypt_exD(num_exD+1:num_exD+8) = bin2_8(1:8); %8bit secret information
                        num_exD = num_exD + 8;
                    else
                        t = num_D - num_exD;
                        Encrypt_exD(num_exD+1:num_exD+t) = bin2_8(1:t); %tbit secret information
                        num_exD = num_exD + t;
                    end
                end
            end
        end
    end
end
end
function [value,this_end] = Huffman_DeCode(Binary,last_end)
% Find the integer value converted into the next Huffman code in the binary bit stream Binary
% Input: Binary (binary mapping sequence), last_end (the position where the previous mapping ends)
% Output: value (decimal integer value)→{0,1,4,5,12,13,14,30,31},end (the end position this time)
len = length(Binary);
i = last_end+1;
t = 0; % count
if i <= len
    if i+1<=len && Binary(i)==0 %→0
        t = t + 1;
        if Binary(i+1) == 0 %→00 means 0
            t = t + 1;
            value = 0;
        elseif Binary(i+1) == 1 %→01 means 1
            t = t + 1;
            value = 1;
        end
    else % Binary(i)==1
        if i+2<=len && Binary(i+1)==0 %→10
            t = t + 2;
            if Binary(i+2) == 0 %→100 means 4
                t = t + 1;
                value = 4;
            elseif Binary(i+2) == 1 %→101 means 5
                t = t + 1;
                value = 5;
            end
        else % Binary(i+1)==1
            if i+3<=len && Binary(i+2)==0 %→110
                t = t + 3;
                if Binary(i+3) == 0 %→1100 means 12
                    t = t + 1;
                    value = 12;
elseif Binary(i+3) == 1 %→1101 means 13
                    t = t + 1;
                    value = 13;
                end
            else % Binary(i+2)==1
                if i+3 <= len
                    t = t + 3;
                    if Binary(i+3) == 0 %→1110 means 14
                        t = t + 1;
                        value = 14;
                    elseif i+4<=len && Binary(i+3)==1 %→1111
                        t = t + 1;
                        if Binary(i+4) == 0 %→11110 means 30
                            t = t + 1;
                            value = 30;
                        elseif Binary(i+4) == 1 %→11111 means 31
                            t = t + 1;
                            value = 31;
                        end
                    else
                        t = 0;
                        value = -1; % indicates that the length of the auxiliary information is not enough to restore the next Huffman code
                    end
                else
                    t = 0;
                    value = -1; % indicates that the length of the auxiliary information is not enough to restore the next Huffman code
                end
            end
        end
    end
else
    t = 0;
    value = -1; % indicates that the length of the auxiliary information is not enough to restore the next Huffman code
end
this_end = last_end + t;
end
function [recover_I] = Recover_Image(stego_I,Image_key,Side_Information,Refer_Value,Map_I,num,ref_x,ref_y)
% Function description: restore the image based on the extracted auxiliary information
% Input: stego_I (secret image), Image_key (image encryption key), Side_Information (auxiliary information), Refer_Value (reference pixel information), Map_I (location map), num (length of secret information), ref_x, ref_y (reference number of rows and columns of pixels)
% output: recover_I (recovery image)
[row,col] = size(stego_I); % count the number of rows and columns of stego_I
%% Restore the reference pixels of the previous ref_y column and the previous ref_x row according to the Refer_Value
refer_I = stego_I;
t = 0; % count
for i=1:row
    for j=1:ref_y
        bin2_8 = Refer_Value(t+1:t+8);
        [value] = Binary_Decimalism(bin2_8); % Convert 8-bit binary array to decimal integer
        refer_I(i,j) = value;
        t = t + 8;
    end
end
for i=1:ref_x
    for j=ref_y+1:col
        bin2_8 = Refer_Value(t+1:t+8);
        [value] = Binary_Decimalism(bin2_8); % Convert 8-bit binary array to decimal integer
        refer_I(i,j) = value;
        t = t + 8;
    end
end
%% Decrypt the image refer_I according to the image encryption key
[decrypt_I] = Encrypt_Image(refer_I,Image_key);
%% Restore pixels at other positions based on Side_Information, Map_I and num
recover_I = decrypt_I;
num_S = length(Side_Information);
num_D = num_S + num; % total number of embedded information 
re = 0; % count
for i=ref_x+1:row
    for j=ref_y+1:col
        if re >= num_D % all bits of embedded information are restored
            break;
        end
        %--------- Find the predicted value of the current pixel---------%
        a = recover_I(i-1,j);
        b = recover_I(i-1,j-1);
        c = recover_I(i,j-1);
        if b <= min(a,c)
            pv = max(a,c);
        elseif b >= max(a,c)
            pv = min(a,c);
        else
            pv = a + c - b;
        end
        %--Convert original and predicted values ​​into 8-bit binary arrays--%
        x = recover_I(i,j);
        [bin2_x] = Decimalism_Binary(x);
        [bin2_pv] = Decimalism_Binary(pv);
        %------- means that this pixel needs to restore 1 bit MSB-------%
        if Map_I(i,j) == 0 %Map=0 means that the 1st MSB of the original pixel value is opposite to its predicted value
            if bin2_pv(1) == 0
                bin2_x(1) = 1;
            else
bin2_x(1) = 0;
            end
            [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
            recover_I(i,j) = value;
            re = re + 1; % restore 1bit
        %------- means that this pixel needs to restore 2 bit MSB-------%
        elseif Map_I(i,j) == 1 %Map=1 means the 2 MSB of the original pixel value is opposite to its predicted value
            if re+2 <= num_D
                if bin2_pv(2) == 0
                    bin2_x(2) = 1;
                else
                    bin2_x(2) = 0;
                end
                bin2_x(1) = bin2_pv(1);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + 2; % restore 2bit
            else
                t = num_D - re; % number of bits remaining to recover
                bin2_x(1:t) = bin2_pv(1:t);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + t; % restore tbit
            end
%------- means that this pixel needs to restore 3 bit MSB-------%
        elseif Map_I(i,j) == 2 %Map=2 means that the 3rd MSB of the original pixel value is opposite to its predicted value
            if re+2 <= num_D
                if bin2_pv(3) == 0
                    bin2_x(3) = 1;
                else
                    bin2_x(3) = 0;
                end
                bin2_x(1:2) = bin2_pv(1:2);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + 3; % restore 3bit
            else
                t = num_D - re; % number of bits remaining to recover
                bin2_x(1:t) = bin2_pv(1:t);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + t; % restore tbit
            end
        %------- means that this pixel needs to restore 4 bit MSB-------%
        elseif Map_I(i,j) == 3 %Map=3 means the 4th MSB of the original pixel value is opposite to its predicted value
            if re+3 <= num_D
                if bin2_pv(4) == 0
                    bin2_x(4) = 1;
                else
                    bin2_x(4) = 0;
                end
bin2_x(1:3) = bin2_pv(1:3);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + 4; % restore 4bit
            else
                t = num_D - re; % number of bits remaining to recover
                bin2_x(1:t) = bin2_pv(1:t);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + t; % restore tbit
            end
        %------- means that this pixel needs to restore 5 bit MSB-------%
        elseif Map_I(i,j) == 4 %Map=4 means that the 5th MSB of the original pixel value is the opposite of its predicted value
            if re+4 <= num_D
                if bin2_pv(5) == 0
                    bin2_x(5) = 1;
                else
                    bin2_x(5) = 0;
                end
                bin2_x(1:4) = bin2_pv(1:4);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + 5; % restore 5bit
            else
                t = num_D - re; % number of bits remaining to recover
                bin2_x(1:t) = bin2_pv(1:t);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + t; % restore tbit
            end
%------- means that this pixel needs to restore 6 bit MSB-------%
        elseif Map_I(i,j) == 5 %Map=5 means the 6MSB of the original pixel value is opposite to its predicted value
            if re+5 <= num_D
                if bin2_pv(6) == 0
                    bin2_x(6) = 1;
                else
                    bin2_x(6) = 0;
                end
                bin2_x(1:5) = bin2_pv(1:5);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + 6; % restore 6bit
            else
                t = num_D - re; % number of bits remaining to recover
                bin2_x(1:t) = bin2_pv(1:t);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + t; % restore tbit
            end
        %------- means that this pixel needs to restore 7 bit MSB-------%
        elseif Map_I(i,j) == 6 %Map=6 means the 7th MSB of the original pixel value is the opposite of its predicted value
            if re+6 <= num_D
                if bin2_pv(7) == 0
                    bin2_x(7) = 1;
                else
                    bin2_x(7) = 0;
                end
                bin2_x(1:6) = bin2_pv(1:6);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + 7; % restore 7bit
            else
t = num_D - re; % number of bits remaining to recover
                bin2_x(1:t) = bin2_pv(1:t);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + t; % restore tbit
            end
        %------- means that this pixel needs to restore 8 bit MSB-------%
        elseif Map_I(i,j) == 7 %Map=7 means the 8th MSB of the original pixel value is opposite to its predicted value
            if re+7 <= num_D
                if bin2_pv(8) == 0
                    bin2_x(8) = 1;
                else
                    bin2_x(8) = 0;
                end
                bin2_x(1:7) = bin2_pv(1:7);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + 8; % restore 8bit
            else
                t = num_D - re; % number of bits remaining to recover
                bin2_x(1:t) = bin2_pv(1:t);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + t; % restore tbit
            end
        %------- means that this pixel needs to restore 8 bit MSB-------%
        elseif Map_I(i,j) == 8 %Map=8 means the original pixel value is equal to its predicted value
            if re+8 <= num_D
                bin2_x(1:8) = bin2_pv(1:8);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + 8; % restore 8bit
            else
t = num_D - re; % number of bits remaining to recover
                bin2_x(1:t) = bin2_pv(1:t);
                [value] = Binary_Decimalism(bin2_x); % Convert 8-bit binary array to decimal integer
                recover_I(i,j) = value;
                re = re + t; % restore tbit
            end
        end
    end
end
end
