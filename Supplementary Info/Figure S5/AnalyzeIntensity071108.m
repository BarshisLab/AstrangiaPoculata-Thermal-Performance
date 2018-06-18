%Coral intensity analysis utility
%Alex Blekhman 2005 (C)

close all
clear all
clc

disp(' ');
disp('*************************************************************************');
disp('*** Coral Color Intensity Analysis Utility 1.0                        ***');
disp('*** Copyright: Alex Blekhman 2005 - 2008 (C)                          ***');
disp('*************************************************************************');
disp(' ');

InputImages=1;
NumPoints=20;
PointWidth=25;
prompt = {'Enter number of input images:','Enter number of points in each image:','Enter averaging window width for each point:'};
dlg_title = 'Input for coral intensity analysis utility';
num_lines= 1;
def     = {num2str(InputImages),num2str(NumPoints),num2str(PointWidth)};
answer  = inputdlg(prompt,dlg_title,num_lines,def);
InputImages=str2num(answer{1});
NumPoints=str2num(answer{2});
PointWidth=str2num(answer{3});
RGB=zeros(NumPoints,3,InputImages);
HSV=zeros(NumPoints,3,InputImages);
for image_num = 1:InputImages
    [filename1, pathname1] = uigetfile( ...
        {'*.jpg;*.jpeg;*.gif;*.bmp','Graphic Files (*.jpg;*.jpeg;*.gif;*.bmp)';
        '*.*',  'All Files (*.*)'}, ...
        ['Pick image file #',num2str(image_num)]);
    if isequal(filename1,0)|isequal(pathname1,0)
        errordlg('File not found','File Error');
    else
        disp(['File ', pathname1, filename1, ' found'])
    end
    im1 = imread([pathname1 filename1]);
    figure('Name',['Image ' int2str(image_num)]);
    imshow(im1);
    title(['Image ' int2str(image_num) ': ' pathname1 filename1]);
    disp(['Please select ',int2str(NumPoints),' points on the coral in image ' int2str(image_num)]);
    % pause
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=1:NumPoints
        figure(image_num)
        [x,y]=ginput(1);
        hold on
        plot(x,y,'g*')
        half_width=(PointWidth-1)/2;
        X=[(x-half_width) (x-half_width) (x+half_width) (x+half_width)];
        Y=[(y-half_width) (y+half_width) (y+half_width) (y-half_width)];
        ribuax=[X, X(1)];  ribuay=[Y, Y(1)];
        plot(ribuax,ribuay,'b-')
        Res1r(1:PointWidth,1:PointWidth)=im1(floor(y-half_width):floor(y+half_width),floor(x-half_width):floor(x+half_width),1);
        Res1_r(i)=median(median(double(Res1r(:,:))));
        Res1g(1:PointWidth,1:PointWidth)=im1(floor(y-half_width):floor(y+half_width),floor(x-half_width):floor(x+half_width),2);
        Res1_g(i)=median(median(double(Res1g(:,:))));
        Res1b(1:PointWidth,1:PointWidth)=im1(floor(y-half_width):floor(y+half_width),floor(x-half_width):floor(x+half_width),3);
        Res1_b(i)=median(median(double(Res1b(:,:))));
    end;
    RGB(:,:,image_num)=[Res1_r;Res1_g;Res1_b]';
    HSV(:,:,image_num)=rgb2hsv(RGB(:,:,image_num));
%     RGB2=[Res2_r;Res2_g;Res2_b]';
%     eval(['RGB',num2str(image_num),'=[Res1_r;Res1_g;Res1_b]']);
%     eval(['HSV',num2str(image_num),'=rgb2hsv(RGB',num2str(i),')']);
end;

%%%%%%%%%%%%%%%RED%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Red intensity profile for selected images');
RED = RGB(:,1,:);
RED = reshape(RED,NumPoints,InputImages);
bar(RED,'group');
legend('Image 1','Image 2');
title('Red  intensity profile analysis for selected images');
xlabel('Point number');
ylabel('Red  intensity coefficient  [0 - darkest, 255 - brightest]');
% disp('*********************************************************************************************************************');
% disp('Red  intensity  full pointwise profile: (row index - number of point, column index - number of picture)');
% disp(RED);

RED_mean=mean(RGB(:,1,:));
RED_mean=reshape(RED_mean,1,InputImages);
figure('Name','RGB RED mean value comparison');
title('RGB RED mean value analysis');
bar(RED_mean,'r');
title('RGB RED mean value comparison');
xlabel('Image number');
ylabel('RGB RED mean value [0 - darkest, 255 - brightest]');
% disp('*********************************************************************************************************************');
% disp('RGB RED mean profile: (column index - number of picture)');
% disp(RED_mean);
%%%%%%%%%%%%%%%RED%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%GREEN%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','GREEN  intensity profile for selected images');
GREEN = RGB(:,2,:);
GREEN = reshape(GREEN,NumPoints,InputImages);
bar(GREEN,'group');
legend('Image 1','Image 2');
title('GREEN  intensity profile analysis for selected images');
xlabel('Point number');
ylabel('GREEN  intensity coefficient  [0 - darkest, 255 - brightest]');
% disp('*********************************************************************************************************************');
% disp('GREEN  intensity  full pointwise profile: (row index - number of point, column index - number of picture)');
% disp(GREEN);

GREEN_mean=mean(RGB(:,2,:));
GREEN_mean=reshape(GREEN_mean,1,InputImages);
figure('Name','RGB GREEN mean value comparison');
title('RGB GREEN mean value analysis');
bar(1:InputImages,GREEN_mean,'g');
title('RGB GREEN mean value comparison');
xlabel('Image number');
ylabel('RGB GREEN mean value [0 - darkest, 255 - brightest]');
% disp('*********************************************************************************************************************');
% disp('RGB GREEN mean profile: (column index - number of picture)');
% disp(GREEN_mean);
%%%%%%%%%%%%%%GREEN%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%BLUE%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','BLUE intensity profile for selected images');
BLUE = RGB(:,3,:);
BLUE = reshape(BLUE,NumPoints,InputImages);
bar(BLUE,'group');
legend('Image 1','Image 2');
title('BLUE  intensity profile analysis for selected images');
xlabel('Point number');
ylabel('BLUE  intensity coefficient  [0 - darkest, 255 - brightest]');
% disp('*********************************************************************************************************************');
% disp('BLUE  intensity  full pointwise profile: (row index - number of point, column index - number of picture)');
% disp(BLUE);

BLUE_mean=mean(RGB(:,3,:));
BLUE_mean=reshape(BLUE_mean,1,InputImages);
figure('Name','RGB BLUE mean value comparison');
title('RGB BLUE mean value analysis');
bar(1:InputImages,BLUE_mean,'b');
title('RGB BLUE mean value comparison');
xlabel('Image number');
ylabel('RGB BLUE mean value [0 - darkest, 255 - brightest]');
% disp('*********************************************************************************************************************');
% disp('RGB BLUE mean profile: (column index - number of picture)');
% disp(BLUE_mean);
%%%%%%%%%%%%%%%BLUE%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%HSV%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure('Name','Full HSV brightness profile for selected images');
% HSV_full = HSV(:,3,:);
% HSV_full = reshape(HSV_full,NumPoints,InputImages);
% bar(HSV_full,'group');
% legend('Image 1','Image 2');
% title('Full brightness profile analysis for selected images');
% xlabel('Point number');
% ylabel('Color brightness coefficient (HSV value) [0 - darkest, 255 - brightest]');
% disp('*********************************************************************************************************************');
% disp('HSV full pointwise profile: (row index - number of point, column index - number of picture)');
% disp(HSV_full);

% HSV_mean=mean(HSV(:,3,:));
% HSV_mean=reshape(HSV_mean,1,InputImages);
% figure('Name','HSV mean value comparison');
% bar(1:InputImages,HSV_mean);
% title('HSV mean value analysis');
% xlabel('Image number');
% ylabel('HSV mean value [0 - darkest, 255 - brightest]');
% disp('*********************************************************************************************************************');
% disp('HSV mean profile: (column index - number of picture)');
% disp(HSV_mean);
%%%%%%%%%%%%%%%HSV%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Final output
% disp('Coral color intensity (Value column - HSV scale):');
% disp('Image 1 Image 2');
% disp(HSV_full);
disp('Coral color intensity (RGB):');
disp('Image 1 [RED, GREEN, BLUE]');
disp([RED(:,1) GREEN(:,1) BLUE(:,1)]);
disp('Image 2 [RED, GREEN, BLUE]');
disp([RED(:,2) GREEN(:,2) BLUE(:,2)]);

disp(' ');
disp('*************************************************************************');
disp('*** Thank you for using Coral Color Intensity Analysis Utility 1.0    ***');
disp('*** Correspondence: Alex Blekhman (ablekhman at gmail dot com)        ***');
disp('*************************************************************************');
disp(' ');
