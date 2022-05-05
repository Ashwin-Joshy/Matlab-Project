[filename, pathname] = uigetfile('*.*', 'Pick an Video');
B=regexp(filename,'\d*','Match')
 disp(B{1})
 filename=strcat(pathname,filename);
