% Coral color calibration utility 1.0
% Copyright by author - Alex Blekhman 2005 - 2007 (C)

function CalibrateImage

close all
clear all
home

warning off all

disp(' ');
disp('******************************************************************');
disp('*** Coral Color Calibration Utility 1.0                        ***');
disp('*** Copyright: Alex Blekhman 2005 - 2008 (C)                   ***');
disp('******************************************************************');
disp(' ');

NumBars = 20;

reply = input('Do you want to use a calibration scale from reference image (R) or to use linear approximation (L)? R/L [L]: ', 's');
if isempty(reply)
    reply = 'L';
end

if(upper(reply) == 'L')
    disp('Performing calibration by using linear scale.');
    full_dynamic_range = 1;
else
    disp('Performing calibration by using reference image.');
    full_dynamic_range = 0;
    reply = input(['Please enter the number of gray level bars on reference image [' int2str(NumBars) ']: '], 's');
    if ~isempty(reply)
        NumBars = str2num(reply);
    end
end

intensity_adjustment = 0;

if(full_dynamic_range == 0)
    disp('Please pick reference image file');
    [filename1, pathname1] = uigetfile( ...
        {'*.jpg;*.jpeg;*.gif;*.bmp','Graphic Files (*.jpg;*.jpeg;*.gif;*.bmp)';
        '*.*',  'All Files (*.*)'}, ...
        'Please pick reference image file');
    if isequal(filename1,0)|isequal(pathname1,0)
        errordlg('File not found','File Error');
        disp('File Error - File not found');
        return;
    else
        disp(['Loading ', pathname1, filename1, ' ...'])
    end
    im1 = imread([pathname1 filename1]);
    disp(['Loaded ', pathname1, filename1, '.'])
    im1_r = im1(:,:,1);im1_g = im1(:,:,2);im1_b = im1(:,:,3);
    im1_r(im1_r == 0) = 1;im1_g(im1_g == 0) = 1;im1_b(im1_b == 0) = 1;
end;

disp('Please pick uncalibrated image file');
[filename2, pathname2] = uigetfile( ...
    {'*.jpg;*.jpeg;*.gif;*.bmp','Graphic Files (*.jpg;*.jpeg;*.gif;*.bmp)';
    '*.*',  'All Files (*.*)'}, ...
    'Please pick uncalibrated image file');
if isequal(filename2,0)|isequal(pathname2,0)
    errordlg('File not found','File Error');
    disp('File Error - File not found');
    return;
else
    disp(['Loading ', pathname2, filename2, ' ...'])
end

im2 = imread([pathname2 filename2]);
disp(['Loaded ', pathname2, filename2, '.'])
im2_r = im2(:,:,1);im2_g = im2(:,:,2);im2_b = im2(:,:,3);
im2_r(im2_r == 0) = 1;im2_g(im2_g == 0) = 1;im2_b(im2_b == 0) = 1;

N1 = 3;
factor = 1.6;

% Get Res1_r, Res1_g, Res1_b - either from a reference image or from a
% linear theoretic scale

if(full_dynamic_range == 1)
    % Create theoretic scale
    Res1_r=[256:-13:0];Res1_g=[256:-13:0];Res1_b=[256:-13:0];
    if(intensity_adjustment == 1)
        Res1_r(N1+1:20) = Res1_r(N1+1:20)/factor;Res1_g(N1+1:20) = Res1_g(N1+1:20)/factor;Res1_b(N1+1:20) = Res1_b(N1+1:20)/factor;
        Res1_r(N1) = Res1_r(N1)/(factor-0.1);Res1_g(N1) = Res1_g(N1)/(factor-0.4);Res1_b(N1) = Res1_b(N1)/(factor-0.4);
        Res1_r(N1-1) = Res1_r(N1-1)/(factor-0.2);Res1_g(N1-1) = Res1_g(N1-1)/(factor-0.5);Res1_b(N1-1) = Res1_b(N1-1)/(factor-.5);
        Res1_r(N1-2) = Res1_r(N1-2)/(factor-0.4);Res1_g(N1-2) = Res1_g(N1-2)/(factor-0.6);Res1_b(N1-2) = Res1_b(N1-2)/(factor-.6);
    end;
else
    figure;
    imshow(im1);
    title(['Reference image: ' pathname1 filename1]);

    width = GetRectangleWidth;

    disp('Please resize reference image so that gray bars are seen clearly and press any key');
    pause
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['Please mark ' int2str(NumBars) ' reference gray bars on the image...']);
    [Res1_r, Res1_g, Res1_b] = GetScaleValues(NumBars,width,im1);
end;

figure;
imshow(im2);
title(['Uncalibrated image: ' pathname2 filename2]);

width = GetRectangleWidth();

disp('Please resize uncalibrated image so that gray bars are seen clearly and press any key');
pause

disp(['Please mark ' int2str(NumBars) ' reference gray bars on the image...']);
[Res2_r, Res2_g, Res2_b] = GetScaleValues(NumBars,width,im2);

if(intensity_adjustment == 0)
    [r1_coeff Res1_r]=LinearRegression(Res1_r');[g1_coeff Res1_g]=LinearRegression(Res1_g');[b1_coeff Res1_b]=LinearRegression(Res1_b');
end;
[r2_coeff Res2_r]=LinearRegression(Res2_r');[g2_coeff Res2_g]=LinearRegression(Res2_g');[b2_coeff Res2_b]=LinearRegression(Res2_b');

r1=[255 Res1_r 0];
g1=[255 Res1_g 0];
b1=[255 Res1_b 0];
r2=[255 Res2_r 0];
g2=[255 Res2_g 0];
b2=[255 Res2_b 0];

figure;
plot( r2 , '-.r');
hold on;
plot( g2 , '-.g');
plot( b2 , '-.b');

plot( r1 , '-ok');
plot( g1 , '-ok');
plot( b1 , '-ok');
title('Calibration parameters - full scale');
legend('Data RED','Data GREEN','Data BLUE','Reference RED, GREEN and BLUE')
xlabel(['Number of gray level strip (1-' int2str(NumBars) ')']);
ylabel('Color intensity (0 - darkest, 255 - brightest)');

% Compute calibrated image
full_r = uint8(round(interp1(r2,r1,[1:256],'linear')));
full_g = uint8(round(interp1(g2,g1,[1:256],'linear')));
full_b = uint8(round(interp1(b2,b1,[1:256],'linear')));

im2_r_new = full_r(im2_r);
im2_g_new = full_g(im2_g);
im2_b_new = full_b(im2_b);

indices = im2_r_new>215&im2_r_new<220&im2_b_new>240;
im2_r_new(indices) = 255;
im2_g_new(indices) = 255;
im2_b_new(indices) = 255;
im2_new = im2;
im2_new(:,:,1) = im2_r_new;
im2_new(:,:,2) = im2_g_new;
im2_new(:,:,3) = im2_b_new;

figure
imshow(im2_new)
title(['Image after calibration from ' filename2])

if(full_dynamic_range == 0)
    figure
    subplot(2,2,1)
    imshow(im1)
    title(['Uncalibrated image: ' filename1]);
    subplot(2,2,2)
    imshow(im2)
    title(['Calibration image: ' filename2])
    subplot(2,2,3)
    imshow(im1)
    title(['Uncalibrated image: ' filename1])
    subplot(2,2,4)
    imshow(im2_new)
    title(['Image after calibration from ' filename2]);
end;

reply = input('Do you want to save the calibrated file? Y/N [Y]: ', 's');
if isempty(reply)
    reply = 'Y';
end

if(upper(reply) == 'Y')
    SaveOutputImage(im2_new);
end;

disp(' ');
disp('******************************************************************');
disp('*** Thank you for using Coral Color Calibration Utility 1.0    ***');
disp('*** Correspondence: Alex Blekhman (ablekhman at gmail dot com) ***');
disp('******************************************************************');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF MAIN FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Res_r, Res_g, Res_b] = GetScaleValues(NumBars,width,im)

for i=1:NumBars,
    %figure(1)
    [x,y]=ginput(1);
    hold on
    half_width = PlotRectangle(x,y,width);
    Resr(1:width,1:width)=im( floor(y-half_width):floor(y+half_width),floor(x-half_width):floor(x+half_width),1);
    Res_r(i)=median(median(double(Resr(:,:))));
    Resg(1:width,1:width)=im( floor(y-half_width):floor(y+half_width),floor(x-half_width):floor(x+half_width),2);
    Res_g(i)=median(median(double(Resg(:,:))));
    Resb(1:width,1:width)=im( floor(y-half_width):floor(y+half_width),floor(x-half_width):floor(x+half_width),3);
    Res_b(i)=median(median(double(Resb(:,:))));
end

function half_width=PlotRectangle(x,y,width)
plot(x,y,'g*')
half_width = (width-1)/2;
X=[(x-half_width) (x-half_width) (x+half_width) (x+half_width)];
Y=[(y-half_width) (y+half_width) (y+half_width) (y-half_width)];
ribuax=[X, X(1)];  ribuay=[Y, Y(1)];
plot(ribuax,ribuay,'b-')

function width = GetRectangleWidth
width = 25;
hold on;
PlotRectangle(width,width,width);
reply = input(['Please enter rectangle width in pixels (see rectangle in upper left corner for example) [' int2str(width) ']: '], 's');
if ~isempty(reply)
    width = str2num(reply);
    if(~isnumeric(width) || width < 0)
        disp('Invalid value for width. Reverting to default (width = 25)');
        width = 25;
    end
end

function SaveOutputImage(im)
disp('Please choose a file name to save the calibrated image');
[filename, pathname, filterindex] = uiputfile( ...
    {'*.jpg;*.jpeg;*.gif;*.bmp','Graphic Files (*.jpg;*.jpeg;*.gif;*.bmp)';
    '*.*',  'All Files (*.*)'}, ...
    'Save as');

if isequal(filename,0) | isequal(pathname,0)
    errordlg('Illegal file','File Error');
    disp('File Error - Illegal file. Could not save the calibrated image');
    return;
end

if(isempty(strfind(filename,'.jpg')) && isempty(strfind(filename,'.jpeg')) && isempty(strfind(filename,'.gif')) && isempty(strfind(filename,'.bmp')))
    disp('File format unrecognized - reverting to default (.jpg)');
    filename = [filename '.jpg'];
end

disp(['Saving the calibrated image to "' filename '".']);

imwrite(im,[pathname filename]);
disp(['Calibrated image saved to ' filename '.']);

function [a,out_p] = LinearRegression(in_p)
size_p = max(size(in_p));
X=[ones(1,size_p)' (1:size_p)'];
a=X\in_p;
out_p=a(1) + a(2)*(1:size_p);
