clear;clc;close all;
I=imread('barcode.jpg');
I2= imread('条形码/IMG_4850.jpg');
gray=rgb2gray(I);
gray2=rgb2gray(I2);
% gausI=imgaussfilt(gray,5);
% bw=double(bw);
% imshow(gausI)
% bw=imbinarize(gausI);
% imshow(bw)
%%
points = detectSURFFeatures(gray); %检测surf特征点
[features, valid_points] = extractFeatures(gray, points); %获取特征点向量(n*64)和特征点坐标

points2=detectSURFFeatures(gray2); %检测surf特征点
[features2, valid_points2] = extractFeatures(gray2, points2); %获取特征点向量(n*64)和特征点坐标
idx = matchFeatures(features,features2);
%%
figure
imshow(gray2)
hold on
plot(valid_points2(idx(:,2)),'showOrientation',true);

location=valid_points2(idx(:,2)).Location; %特征点的位置
s1 = std2(location(:,1)); %计算均方差
s2=std2(location(:,2)); 
%%
gausI=imgaussfilt(gray2,15); %高斯过滤
bw=gausI>120 & gausI<150; %阈值筛选
% bw2=conv2(bw,ones(5),'same');
% bw2=bw2>.5;
figure
imshow(bw)

%%
% b2=zeros(size(bw));
% b2(clist{1})=1;
% imshow(b2)
%%
conn=bwconncomp(bw); %找联通区域
clist=conn.PixelIdxList;
maxC=0;
for i=1:length(clist)
    if numel(clist{i})>1000
        [m,n]=ind2sub(size(bw),clist{i}); %下标转换为坐标
        if std(m)<s2*2 && std(n)<s1*2
            %在当前连通域里是否有surf特征点
            mpos=find(location(:,1)>min(n) & location(:,1)<max(n)); 
            npos=find(location(:,2)>min(m) & location(:,2)<max(m));
            c=intersect(mpos,npos);
            if numel(c)>maxC
                maxC=numel(c);
                labelid=i;
                final_m=m;
                final_n=n;
            end
        end
    end
end

newI=I2(min(final_m):max(final_m),min(final_n):max(final_n),:);
figure
imshow(newI)
%%
bw2=rgb2gray(newI)<60;
bw3=bwareaopen(bw2,150);
imshow(bw3)
%%
minL=99999;
for i=-30:.5:30
    tempI=imrotate(bw3,i,'loose');
    [m,n]=find(tempI);
    L=length(unique(n))+length(unique(m));
    if L<minL
        minL=L;
        ang=i;
        tempBW=tempI;
    end
end
%%
bw4=conv2(tempBW,ones(1,27),'same');
conn=bwconncomp(bw4);
clist=conn.PixelIdxList;
L=0;
for i=1:length(clist)
    if numel(clist{i})>10000
        [m,n]=ind2sub(size(bw4),clist{i});
        mT=max(m)^2+max(n)^2; %找最接近右下角的连通域
        if mT>L
            L=mT;
            cidx=i;
            f_m=m;
            f_n=n;
        end
    end
end
%%

newI2=imrotate(newI,ang,'loose');
newI3=newI2(min(f_m):max(f_m),min(f_n):max(f_n),:);
figure
% J=imadjust(newI3); %增加对比度
imshow(newI3)
ocrResults = ocr(newI3);
title(ocrResults.Words)