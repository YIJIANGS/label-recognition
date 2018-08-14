clear;clc;close all;
[file,path,indx] = uigetfile('*');
im=imread([path,file]); %读取图片
figure
imshow(im)
if numel(im)<1e7 %判断图像是否是小图（筛选不同阈值）
    isSmall=1;
else
    isSmall=0;
end
hold on
title('原图')
I2=imread('label.jpg');
points = detectSURFFeatures(rgb2gray(I2)); %检测surf特征点
[features, valid_points] = extractFeatures(rgb2gray(I2), points); %获取特征点向量(n*64)和特征点坐标
points2 = detectSURFFeatures(rgb2gray(im)); %检测surf特征点
[features2, valid_points2] = extractFeatures(rgb2gray(im), points2); %获取特征点向量(n*64)和特征点坐标
idx = matchFeatures(features,features2);
plot(valid_points2(idx(:,2)),'showOrientation',true);
location=valid_points2(idx(:,2)).Location; %特征点的位置
if numel(location)>4 %确定有标签
%分别提取rgb三色通道
r=im(:,:,1); 
g=im(:,:,2);
b=im(:,:,3);
%找出蓝色标签范围（手动测试）
r1=r<50;
g1=g>50 & g<160;
b1=b>70;
bw=r1.*g1.*b1; %同时符合三色范围的点    
    
%% 提取大块标识
bw2=conv2(bw,ones(2),'same'); %卷积扩充标识标识范围，作用相当于膨胀
bw3=imfill(bw2,'holes'); %填充膨胀图
conn=bwconncomp(bw3); %找到联通区域
plist=conn.PixelIdxList; %所有联通区域下标
pnum=cellfun(@length,plist); %所有联通区域内联通点数量
[~,idx]=max(pnum); %找出最大联通区域所在
%重新构建三色通道 所有非标识区域涂黑
r3=zeros(size(r1)); g3=zeros(size(r1)); b3=zeros(size(r1));
% r3=nan(size(r1)); g3=nan(size(r1)); b3=nan(size(r1));
%所有单色通道重新上色
r3(plist{idx})=r(plist{idx});
g3(plist{idx})=g(plist{idx});
b3(plist{idx})=b(plist{idx});
im2=uint8(cat(3,r3,g3,b3)); %三色通道融合
[m,n]=find(rgb2gray(im2)); %找到蓝色区域坐标
im2=im(min(m):max(m),min(n):max(n),:); %截取蓝色区域图
% figure
% imshow(im2) %显示大块标识图像
% title('大标识图')
%% 大标识旋转矫正
bw4=bw3(min(m):max(m),min(n):max(n),:); %从原矩阵二值图中截取大标识区域
% bw4=zeros(size(bw3));
% bw4(plist{idx})=1;
bw4 = bwconvhull(bw4); %通过边缘点拟合多边形
minL=99999;
%设置一个角度范围 以每0.5的角度进行图像旋转，每次旋转后统计其水平和垂直方向的像素点
%返回两个方向占用最少的角度
for i=-30:.5:30
    tempI=imrotate(bw4,i,'loose');
    [m,n]=find(tempI);
    L=max(n)-min(n)+max(m)-min(m);
    if L<minL
        minL=L;
        ang=i;
    end
end
im3=imrotate(im2,ang,'crop'); %确定角度后对图像进行旋转
% imshow(im3)
%% 提取小块标识
%分别提取rgb三色通道
r=im3(:,:,1); 
g=im3(:,:,2);
b=im3(:,:,3);
%找出白色标签范围（手动测试）
r1=r>180;
g1=g>170;
b1=b>150;
bw5=r1.*g1.*b1; %同时符合三色范围的点
bw6=conv2(bw5,ones(2),'same'); %卷积扩充标识标识范围，作用相当于膨胀
bw7=imfill(bw6,'holes'); %填充空隙
conn=bwconncomp(bw7); %找到联通区域
plist=conn.PixelIdxList; %所有联通区域下标
pnum=cellfun(@length,plist); %所有联通区域内联通点数量
[~,idx]=sort(pnum,'descend'); %找出最大2个联通区域所在

%分别把2块最大的联通区域的下标转换为坐标
[x,y]=ind2sub(size(bw7),plist{idx(1)}); 
[x1,y1]=ind2sub(size(bw7),plist{idx(2)});
figure
imshow(im3)
hold on
%分别找到2块区域的边界
k = boundary(x,y,.02); %拟合边界
plot(y(k),x(k),'r') 
k=boundary(x1,y1,.02);
plot(y1(k),x1(k),'r')
if mean(x)<mean(x1) %通过y值平均值判断标签的位置
    idx1=idx(1);
    idx2=idx(2);
else
    idx1=idx(2);
    idx2=idx(1);
end
title('大标识旋转后与小标识定位')
%% 小标识区域
bw8=zeros(size(bw7)); %创建一个全黑矩阵
bw8(plist{idx1})=1;
bw9 = bwconvhull(bw8); %边界拟合
[m,n]=find(bw9);
%重新构建三色通道 所有非标识区域涂黑
r3=zeros(size(r1)); g3=zeros(size(r1)); b3=zeros(size(r1));
%所有单色通道重新上色
r3(bw9)=r(bw9);
g3(bw9)=g(bw9);
b3(bw9)=b(bw9);
im4=uint8(cat(3,r3,g3,b3)); %三色通道融合 并转为uint8形
im4=im3(min(m):max(m),min(n):max(n),:);
% imshow(im4) %显示小块标识图像
%% 小标识旋转矫正
minL=1000;
for i=-30:.5:30
    tempI=imrotate(bw9,i,'loose');
    [m,n]=find(tempI);
    L=length(unique(n))+length(unique(m));
    if L<minL
        minL=L;
        ang=i;
    end
end
im5=imrotate(im4,ang,'crop');

figure
imshow(im5)

%% 能效图区域
bw8=zeros(size(bw7));
bw8(plist{idx2})=1;
bw9 = bwconvhull(bw8);
[m,n]=find(bw9);
%重新构建三色通道 所有非标识区域涂黑
% r3=zeros(size(r)); g3=zeros(size(r)); b3=zeros(size(r));
% %所有单色通道重新上色
% r3(bw9)=r(bw9);
% g3(bw9)=g(bw9);
% b3(bw9)=b(bw9);
% im4=uint8(cat(3,r3,g3,b3)); %三色通道融合 并转为uint8形
im4=im3(min(m):max(m),min(n):max(n),:);
minL=99999;
for i=-30:.5:30
    tempI=imrotate(bw9,i,'loose');
    [m,n]=find(tempI);
    L=length(unique(n))+length(unique(m));
    if L<minL
        minL=L;
        ang=i;
    end
end
im6=imrotate(im4,ang,'crop');
figure
imshow(im6)
title('能效图')
%% 识别能效图

gray=rgb2gray(im6);
bw10=gray<40 & gray>0;
bw11=conv2(bw10,ones(1,30),'same');
bw11=bwareaopen(bw11,300);
% imshow(bw11)
% PSF = fspecial('gaussian',11,2.1);
% J1 = deconvlucy(gray,PSF);
% ocrResults=ocr(J1);
% ocrResults.Words
% imshow(J1)
%%
conn=bwconncomp(bw11);
clist=conn.PixelIdxList;
L=0;
for i=1:length(clist)
    [m,n]=ind2sub(size(bw11),clist{i});
    mT=max(n)+(size(bw11,1)-min(m)); %找最接近右上角的连通域
    if mT>L
        L=mT;
        cidx=i;
        f_m=m;
        f_n=n;
    end
end
im7=im6(min(f_m):max(f_m),min(f_n):max(f_n),:);
%%
ocrResults=ocr(im6);
wpos=ocrResults.WordBoundingBoxes;
[~,idx]=min(abs(min(f_m)-wpos(:,2))+abs(min(f_n)-wpos(:,1)));
figure
imshow(im7)
if ~isempty(ocrResults.Words)
title(ocrResults.Words{idx})
end
fileID = fopen('result.txt','w');
fprintf(fileID,['识别的结果为： ',ocrResults.Words{idx}]);
fclose(fileID);
else
    title('没有能效标签')
end
% title(ocrResults.Words)
%% 二维码
qr=imread('qrcode.jpg');
qr=rgb2gray(qr);
points = detectSURFFeatures(qr); %检测surf特征点
[features, valid_points] = extractFeatures(qr, points); %获取特征点向量(n*64)和特征点坐标
idx = matchFeatures(features,features2);
figure
imshow(im)
hold on
plot(valid_points2(idx(:,2)),'showOrientation',true);
title('二维码检测')
location=valid_points2(idx(:,2)).Location; %特征点的位置
if numel(location)>10
s1 = std2(location(:,1)); %计算均方差
s2=std2(location(:,2)); 
%%
% if isSmall==1 %如果是小图 设置小滤窗口
% gausI=imgaussfilt(rgb2gray(im),5); %高斯过滤
% else
%     gausI=imgaussfilt(rgb2gray(im),15);
% end
bw=rgb2gray(im)<20; %阈值筛选
if isSmall==1
    bw2=bwmorph(bw,'dilate',15);
else
bw2=bwmorph(bw,'dilate',25); %膨胀处理
end
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
figure
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
figure
imshow(newI2)
else
    title('未检测到二维码')
end