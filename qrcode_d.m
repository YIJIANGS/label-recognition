clear;clc;close all;
qr=imread('qrcode.jpg');
qr=rgb2gray(qr);
points = detectSURFFeatures(qr); %检测surf特征点
[features, valid_points] = extractFeatures(qr, points); %获取特征点向量(n*64)和特征点坐标
%%
im=imread('二维码2/IMG_4891.jpg');
gray=rgb2gray(im);
points2 = detectSURFFeatures(gray); %检测surf特征点
[features2, valid_points2] = extractFeatures(gray, points2); %获取特征点向量(n*64)和特征点坐标
idx = matchFeatures(features,features2);
%%
imshow(im)
hold on
plot(valid_points2(idx(:,2)),'showOrientation',true);
location=valid_points2(idx(:,2)).Location; %特征点的位置
s1 = std2(location(:,1)); %计算均方差
s2=std2(location(:,2)); 
%%
gausI=imgaussfilt(gray,15); %高斯过滤
bw=gausI<60; %阈值筛选
bw2=bwmorph(bw,'dilate',35); %膨胀处理
figure
imshow(bw2)

%%
conn=bwconncomp(bw2); %找联通区域
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

newI=im(min(final_m):max(final_m),min(final_n):max(final_n),:);
figure
imshow(newI)
%%
bw2=rgb2gray(newI)<40;
bw3=bwareaopen(bw2,150); %去噪处理
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
newI2=imrotate(newI,ang,'loose');
imshow(newI2)