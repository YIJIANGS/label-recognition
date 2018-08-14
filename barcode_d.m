clear;clc;close all;
I=imread('barcode.jpg');
I2= imread('������/IMG_4850.jpg');
gray=rgb2gray(I);
gray2=rgb2gray(I2);
% gausI=imgaussfilt(gray,5);
% bw=double(bw);
% imshow(gausI)
% bw=imbinarize(gausI);
% imshow(bw)
%%
points = detectSURFFeatures(gray); %���surf������
[features, valid_points] = extractFeatures(gray, points); %��ȡ����������(n*64)������������

points2=detectSURFFeatures(gray2); %���surf������
[features2, valid_points2] = extractFeatures(gray2, points2); %��ȡ����������(n*64)������������
idx = matchFeatures(features,features2);
%%
figure
imshow(gray2)
hold on
plot(valid_points2(idx(:,2)),'showOrientation',true);

location=valid_points2(idx(:,2)).Location; %�������λ��
s1 = std2(location(:,1)); %���������
s2=std2(location(:,2)); 
%%
gausI=imgaussfilt(gray2,15); %��˹����
bw=gausI>120 & gausI<150; %��ֵɸѡ
% bw2=conv2(bw,ones(5),'same');
% bw2=bw2>.5;
figure
imshow(bw)

%%
% b2=zeros(size(bw));
% b2(clist{1})=1;
% imshow(b2)
%%
conn=bwconncomp(bw); %����ͨ����
clist=conn.PixelIdxList;
maxC=0;
for i=1:length(clist)
    if numel(clist{i})>1000
        [m,n]=ind2sub(size(bw),clist{i}); %�±�ת��Ϊ����
        if std(m)<s2*2 && std(n)<s1*2
            %�ڵ�ǰ��ͨ�����Ƿ���surf������
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
        mT=max(m)^2+max(n)^2; %����ӽ����½ǵ���ͨ��
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
% J=imadjust(newI3); %���ӶԱȶ�
imshow(newI3)
ocrResults = ocr(newI3);
title(ocrResults.Words)