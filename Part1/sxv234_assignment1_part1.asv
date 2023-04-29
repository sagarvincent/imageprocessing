%take in the image 
malard = imread('Malards.jpg');

%seperate the red channel 
red_ch = malard(:,:,1);

load roberts.mat

x = conv2(robertsA,red_ch);
y = conv2(robertsB,red_ch);

mag = x.*x + y.*y;
new = sqrt(mag);
show_image(abs(new)>50)
title('Generating edges > 50 using roberts filter')

%%
%first filter
mat1 = [[0,-1,-1,-1,-1,0]
    [-1,-2,-3,-3,-2,-1]
    [-1,-3,12,12,-2,-1]
    [-1,-3,12,12,-3,-1]
    [-1,2,-3,-3,-2,-1]
    [0,-1,-1,-1,-1,0]];
%second filter
mat2 = [[1,1,1,1,1,1]
    [-1,-2,-3,-3,-2,-1]
    [-1,-3,-4,-4,-3,-1]
    [1,3,4,4,3,1]
    [1,2,3,3,2,1]
    [-1,-1,-1,-1,-1,-1]];

mat1,mat2

m = rgb2gray(malard);

%convolve using the first filter
c1 = convn(mat1,m);
%convolve using the second filter
c2 = convn(mat2,m);
%display the results
show_image(c1)
%title('Q 1.2 : First Filter image')
show_image(c2)
%title('Q 1.2 : Second image')\

%convolve using the first filter
c1T = convn(mat1.',m);
%convolve using the second filter
c2T = convn(mat2.',m);
%display the results
show_image(c1T)
%title('Q 1.2 : First Filter image')
show_image(c2T)
%title('Q 1.2 : Second image')
%%
earth = imread(['D:\Practical\Matlab_ws\Robot Vision [06-25024] Summative_1\Part 1\','Earth.jpg']);
fil1 = DiffGaus(im2gray(earth),15,1,2*sqrt(2),1);
fil = DiffGaus(earth,15,1,2*sqrt(2),0);
show_image(earth)
show_image(fil)
title('')
show_image(fil1)
title('Laplacian')
%%
I1 = earth
im = uint8.empty(0,5);
k=1
figure();
for i = 1:4    
    T1 = im2gray(I1) ;
    im{i} = imresize(T1,k);
    a = cell2mat(im(i));
    k = k/2;
    I1 = T1;
end

im2 = uint8.empty(0,5);

for i = 1:4
    I2 = cell2mat(im(i));
    t = DiffGaus(I2,15,1,2*sqrt(2),0);
    im2{i} = t;    
end
for i = 1:4
    a = cell2mat(im2(i));
    show_image(a);
    title('Guassian without binary')

end

%% 
%
fig = ones(1456,751,3);
r = fig(:,:,1);
g = fig(:,:,2) ;
b = fig(:,:,3);
[x,y] = size(r);
for i = 1:4
    a = cell2mat(im2(i));
    %converting to binary
    ra = im2bw(a,.01); 
    [xa,ya] = size(ra);
    add = zeros(376,1);
    
    if (y-ya) ==  0
        r(x-xa:x,y-ya:y) = ra;
        x = x -xa;

    else
        y1 = floor((y-ya)/2);
        y2 = ceil((y-ya)/2);
        if y1 == 0
            y1 = 1;
            y2 = 0;
        else
            y1 = y2;
        end
        
        if ya == 375
            ra1 = [ra add];
        else
            ra1 = ra;
        end
        
        if x ==xa
            xa = 47;
            x = x+2;
            y= y -1;

        end
        r(x-xa:x-1,y1:y-y1) = ra1;
        x = x - xa;
        
       

    end
end
f = r; 
show_image(f)
%%

%%

%%
function y = DiffGaus(I,fsize,sig1,sig2,lap)
    n = fsize;
    f1 = zeros(fsize,fsize);
    for i = 1:fsize
        for j = 1:fsize
            f1(i,j) = 1/(2*pi*sig1*sig1) * exp(-((i-(n+1)/2)^2+(j-(n+1)/2)^2)/(2*sig1^2));
        end
    end
    f2 = zeros(fsize,fsize);
    for i = 1:fsize
        for j = 1:fsize
            f2(i,j) = 1/(2*pi*sig2*sig2) * exp(-((i-(n+1)/2)^2+(j-(n+1)/2)^2)/(2*sig2^2));
        end
    end
    
    y1 = imfilter(I,f1);
    y2 = imfilter(y1,f2);
    y2 = y2 - y1;

    y2 = uint8(y2);
    if lap==1
        y2 = I;
        lap = fspecial('laplacian');
        y = imfilter(y1,lap);
        %y = locallapfilt(y2,2,.2);
    else
        y = im2gray(y2); 
    end

end