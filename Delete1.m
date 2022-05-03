app=1
%Convert_To_Video(app)
%Convert_To_Frames(app)  
decryption()
function []=Convert_To_Video(app)
           %This is a function to convert encrypted images to video 
           outputVideo = VideoWriter('Video23.avi','Uncompressed AVI');
           outputVideo.FrameRate = 30;%Setting the frame rate  
           open(outputVideo)
           encryptedImageFolder='FinalImages'
           encryptedImagePath=dir([encryptedImageFolder '/*.tiff']);
           numberOfEncryptedImages=size(encryptedImagePath,1);
           %Looping through all images and adding to video
           for i = 1:numberOfEncryptedImages
            img = imread(['FinalImages\finalImage' num2str(i) '.tiff']);
            writeVideo(outputVideo,uint8(img))
           end
           close(outputVideo)
        end
function []=Convert_To_Frames(app)
            %Step for choosing the video
            [filename, pathname] = uigetfile('D:\project\Main Project\Matlab-Project/*.*', 'Pick an Image');
            filename=strcat(pathname,filename);
            %Changing current folder
            cd (pathname)
            cd ..
            EncryptedVideo=VideoReader(filename);
            vid = read(EncryptedVideo);
             % read the total number of frames
            frames = EncryptedVideo.NumFrames;
            % file format of the frames to be saved in
            ST ='.tiff';
            % reading and writing the frames 
             k=[]
             for x = 1 : 90%frames
              
                % converting integer to string
                Sx = num2str(x);
                % concatenating 2 strings
                Strc = strcat(Sx, ST);
                Vid = vid(:, :, :, x);
                cd frames\encryptedVideoFrames
               
                % exporting the frames
                pause(0.000001);
                imwrite(Vid, Strc);
                cd ..\..
             end
end  
function [] = decryption()
             Mosaic=imread('DecryptedImages\mosaic1.tiff');
             n=4;
            [M,N,ch]=size(Mosaic);
            numberOfblocks=M*N/n^2;
            Mosaic=double(Mosaic(:));
            len_ld=Mosaic(end);
            for i=len_ld:-1:1
               ld(i)=num2str(Mosaic(end-i));
            end
            ld=str2double(ld);
            dat=rem(Mosaic,4);
            dat=dat(1:ld);
            dat=de2bi(dat,2);
            data=dat(:);
        dec=bi2de(reshape(data,length(data)/8,8))';
        muR1=dec(1:numberOfblocks);
        muG1=dec(numberOfblocks+1:2*numberOfblocks);
        muB1=dec(2*numberOfblocks+1:3*numberOfblocks);
        muR=dec(3*numberOfblocks+1:4*numberOfblocks);
        muG=dec(4*numberOfblocks+1:5*numberOfblocks);
        muB=dec(5*numberOfblocks+1:6*numberOfblocks);
        avgstd=dec(6*numberOfblocks+1:7*numberOfblocks);
        avgstd1=dec(7*numberOfblocks+1:8*numberOfblocks);
        [val1,tileindx]=sort(avgstd);
        [val2,targetindx]=sort(avgstd1);
        Mosaic=reshape(Mosaic,[M,N,3]);
        k=1;
        for i=1:n:M
            for j=1:n:N
                ntile{k}=Mosaic(i:i+n-1,j:j+n-1,:);
                k=k+1;
            end
        end
        %secret Image Extraction
        ntile=ntile(targetindx);
        for k=1:numberOfblocks
            ci=ntile{k};

           nci(:,:,1)=(ci(:,:,1)-muR1(k))+muR(k);
           nci(:,:,2)=(ci(:,:,2)-muG1(k))+muG(k);
           nci(:,:,3)=(ci(:,:,3)-muB1(k))+muB(k);
           ntile3{k}=nci;
        end
        ntile3(tileindx)=ntile3;
        k=1;
        for i=1:n:M
            for j=1:n:N
               RI(i:i+n-1,j:j+n-1,:)=ntile3{k};
                k=k+1;
            end
        end
        figure,imshow(uint8(RI));
        title('Retrieved Secret')
        avg=0;
end       