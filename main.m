clear;clc;close all;
[file,path,indx] = uigetfile('*');
im=imread([path,file]); %读取图片
% im=imread('能识别图片/IMG_4918(1).jpg');
im=imresize(im,[2000 nan]);
tic
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
%%
if numel(location)>2 %确定有标签
%分别提取rgb三色通道
r=im(:,:,1); 
g=im(:,:,2);
b=im(:,:,3);
%找出蓝色标签范围（手动测试）
r1=r<50;
g1=g>50 & g<160;
b1=b>70;
bw=r1.*g1.*b1; %同时符合三色范围的点    
    
%% 提取大块标识并矫正
[ul,ur,lr,ll]=find_4_corners(bw);
im3=perspective_transform(im,bw,ul,ur,lr,ll);

% imshow(im3)
% 提取小块标识
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
% k = boundary(x,y,.02); %找到边界点
% plot(y(k),x(k),'r') %下标转坐标的过程中有倒置，所以yx
% k=boundary(x1,y1,.02);
% plot(y1(k),x1(k),'r')
if mean(x)<mean(x1) %通过y值平均值判断标签的位置
    idx1=idx(1);%idx（）代表在px里的位置
    idx2=idx(2);
else
    idx1=idx(2);
    idx2=idx(1);
end
title('大标识矫正')
%% 小标识区域
bw8=zeros(size(bw7)); %创建一个全黑矩阵
bw8(plist{idx1})=1;
bw9 = bwconvhull(bw8); %边界拟合
[m,n]=find(bw9);
%重新构建三色通道 所有非标识区域涂黑
% %r3=zeros(size(r1)); g3=zeros(size(r1)); b3=zeros(size(r1));
% %所有单色通道重新上色
% %r3(bw9)=r(bw9);
% g3(bw9)=g(bw9);
% b3(bw9)=b(bw9);
% im4=uint8(cat(3,r3,g3,b3)); %三色通道融合 并转为uint8形
im4=im3(min(m):max(m),min(n):max(n),:);
% imshow(im4) %显示小块标识图像

%%
g=rgb2gray(im4);
g=imadjust(g);
% imshow(g)
b2=g>220;
b3=bwareaopen(b2,1000);
b3=bwconvhull(b3);
b4=imcomplement(b2).*b3;
figure
imshow(b4)
hold on
conn=bwconncomp(b4);
clist3=conn.PixelIdxList;
pnum=cellfun(@length,clist3);
[pn,idx]=sort(pnum,'descend');
if pn(5)+pn(6)<pn(4)
    for i=1:4
        [m,n]=ind2sub(size(b2),clist3{idx(i)});
        meanM(i)=mean(m);
        rectangle('position',[min(n),min(m),max(n)-min(n),max(m)-min(m)],'edgecolor','r')
    end
    d=pdist(meanM');
        dm=[1 1 1 2 2 3; 2 3 4 3 4 4];
        Midx=dm(:,find(d==min(d)));
        meanM(Midx(2))=[];
        [newOrder,orderId]=sort(meanM);
        label=['能耗有3级，检测出来的级数为',num2str(find(orderId==Midx(1)))];
else
   for i=1:6
        [m,n]=ind2sub(size(b2),clist3{idx(i)});
        meanM(i)=mean(m);
        rectangle('position',[min(n),min(m),max(n)-min(n),max(m)-min(m)],'edgecolor','r')
    end
    d=pdist(meanM');
        dm=[1 1 1 1 1 2 2 2 2 3 3 3 4 4 5; 2 3 4 5 6 3 4 5 6 4 5 6 5 6 6];
        Midx=dm(:,find(d==min(d)));
        meanM(Midx(2))=[];
        [newOrder,orderId]=sort(meanM);
        label=['能耗有5级，检测出来的级数为',num2str(find(orderId==Midx(1)))];
end
title(label)
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
im6=im3(min(m):max(m),min(n):max(n),:);

figure
imshow(im6)
title('能效图')

%% 识别能效图

gray=rgb2gray(im6);
bw10=gray<40 & gray>0;%二值化
if isSmall==1
bw11=conv2(bw10,ones(1,30),'same');
else
bw11=conv2(bw10,ones(1,75),'same');  
end
bw11=bwareaopen(bw11,200);%区噪声
% imshow(bw11)

conn=bwconncomp(bw11);%找到白色连通域
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

ocrResults=ocr(im6);%ocr识别im6所有字符
wpos=ocrResults.WordBoundingBoxes;%找到识别的位置
[~,idx]=min(abs(min(f_m)-wpos(:,2))+abs(min(f_n)-wpos(:,1)));%找到ocr里最接近的位置
figure
imshow(im7)
try
title(ocrResults.Words{idx})
fileID = fopen('result.txt','w');
fprintf(fileID,['识别的结果为： ',ocrResults.Words{idx}]);
fclose(fileID);
catch
end
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
if numel(location)>2
s1 = std2(location(:,1)); %计算xy均方差
s2=std2(location(:,2)); 
%%
% if isSmall==1 %如果是小图 设置小滤窗口
% gausI=imgaussfilt(rgb2gray(im),5); %高斯过滤,找到黑色部分
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
%         if std(m)<s2*2 && std(n)<s1*2
            %在当前连通域里是否有surf特征点
            mpos=find(location(:,1)>min(n) & location(:,1)<max(n)); 
            npos=find(location(:,2)>min(m) & location(:,2)<max(m));
            c=intersect(mpos,npos);%找到同时符合mn方向的特征点
            if numel(c)>maxC%统计特征点最多的区域
                maxC=numel(c);
                labelid=i;
                final_m=m;
                final_n=n;
            end
%         end
    end
end

newI=im(min(final_m):max(final_m),min(final_n):max(final_n),:);
figure
imshow(newI)

%%
bw4=zeros(size(bw2));
bw4(min(final_m):max(final_m),min(final_n):max(final_n))=bw2(min(final_m):max(final_m),min(final_n):max(final_n));
bw4=bwareaopen(bw4,1000);
bw4=double(bw4);
[ul,ur,lr,ll]=find_4_corners(bw4);
newI2=perspective_transform(im,bw4,ul,ur,lr,ll);
figure
imshow(newI2)
else
    title('未检测到二维码')
end
toc