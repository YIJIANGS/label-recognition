clear;clc;close all;
[file,path,indx] = uigetfile('*');
im=imread([path,file]); %��ȡͼƬ
% im=imread('��ʶ��ͼƬ/IMG_4918(1).jpg');
im=imresize(im,[2000 nan]);
tic
figure
imshow(im)
if numel(im)<1e7 %�ж�ͼ���Ƿ���Сͼ��ɸѡ��ͬ��ֵ��
    isSmall=1;
else
    isSmall=0;
end
hold on
title('ԭͼ')
I2=imread('label.jpg');
points = detectSURFFeatures(rgb2gray(I2)); %���surf������
[features, valid_points] = extractFeatures(rgb2gray(I2), points); %��ȡ����������(n*64)������������
points2 = detectSURFFeatures(rgb2gray(im)); %���surf������
[features2, valid_points2] = extractFeatures(rgb2gray(im), points2); %��ȡ����������(n*64)������������
idx = matchFeatures(features,features2);
plot(valid_points2(idx(:,2)),'showOrientation',true);
location=valid_points2(idx(:,2)).Location; %�������λ��
%%
if numel(location)>2 %ȷ���б�ǩ
%�ֱ���ȡrgb��ɫͨ��
r=im(:,:,1); 
g=im(:,:,2);
b=im(:,:,3);
%�ҳ���ɫ��ǩ��Χ���ֶ����ԣ�
r1=r<50;
g1=g>50 & g<160;
b1=b>70;
bw=r1.*g1.*b1; %ͬʱ������ɫ��Χ�ĵ�    
    
%% ��ȡ����ʶ������
[ul,ur,lr,ll]=find_4_corners(bw);
im3=perspective_transform(im,bw,ul,ur,lr,ll);

% imshow(im3)
% ��ȡС���ʶ
%�ֱ���ȡrgb��ɫͨ��
r=im3(:,:,1); 
g=im3(:,:,2);
b=im3(:,:,3);
%�ҳ���ɫ��ǩ��Χ���ֶ����ԣ�
r1=r>180;
g1=g>170;
b1=b>150;
bw5=r1.*g1.*b1; %ͬʱ������ɫ��Χ�ĵ�
bw6=conv2(bw5,ones(2),'same'); %���������ʶ��ʶ��Χ�������൱������
bw7=imfill(bw6,'holes'); %����϶
conn=bwconncomp(bw7); %�ҵ���ͨ����
plist=conn.PixelIdxList; %������ͨ�����±�
pnum=cellfun(@length,plist); %������ͨ��������ͨ������
[~,idx]=sort(pnum,'descend'); %�ҳ����2����ͨ��������

%�ֱ��2��������ͨ������±�ת��Ϊ����
[x,y]=ind2sub(size(bw7),plist{idx(1)}); 
[x1,y1]=ind2sub(size(bw7),plist{idx(2)});
figure
imshow(im3)
hold on
%�ֱ��ҵ�2������ı߽�
% k = boundary(x,y,.02); %�ҵ��߽��
% plot(y(k),x(k),'r') %�±�ת����Ĺ������е��ã�����yx
% k=boundary(x1,y1,.02);
% plot(y1(k),x1(k),'r')
if mean(x)<mean(x1) %ͨ��yֵƽ��ֵ�жϱ�ǩ��λ��
    idx1=idx(1);%idx����������px���λ��
    idx2=idx(2);
else
    idx1=idx(2);
    idx2=idx(1);
end
title('���ʶ����')
%% С��ʶ����
bw8=zeros(size(bw7)); %����һ��ȫ�ھ���
bw8(plist{idx1})=1;
bw9 = bwconvhull(bw8); %�߽����
[m,n]=find(bw9);
%���¹�����ɫͨ�� ���зǱ�ʶ����Ϳ��
% %r3=zeros(size(r1)); g3=zeros(size(r1)); b3=zeros(size(r1));
% %���е�ɫͨ��������ɫ
% %r3(bw9)=r(bw9);
% g3(bw9)=g(bw9);
% b3(bw9)=b(bw9);
% im4=uint8(cat(3,r3,g3,b3)); %��ɫͨ���ں� ��תΪuint8��
im4=im3(min(m):max(m),min(n):max(n),:);
% imshow(im4) %��ʾС���ʶͼ��

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
        label=['�ܺ���3�����������ļ���Ϊ',num2str(find(orderId==Midx(1)))];
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
        label=['�ܺ���5�����������ļ���Ϊ',num2str(find(orderId==Midx(1)))];
end
title(label)
%% ��Чͼ����
bw8=zeros(size(bw7));
bw8(plist{idx2})=1;
bw9 = bwconvhull(bw8);
[m,n]=find(bw9);
%���¹�����ɫͨ�� ���зǱ�ʶ����Ϳ��
% r3=zeros(size(r)); g3=zeros(size(r)); b3=zeros(size(r));
% %���е�ɫͨ��������ɫ
% r3(bw9)=r(bw9);
% g3(bw9)=g(bw9);
% b3(bw9)=b(bw9);
% im4=uint8(cat(3,r3,g3,b3)); %��ɫͨ���ں� ��תΪuint8��
im6=im3(min(m):max(m),min(n):max(n),:);

figure
imshow(im6)
title('��Чͼ')

%% ʶ����Чͼ

gray=rgb2gray(im6);
bw10=gray<40 & gray>0;%��ֵ��
if isSmall==1
bw11=conv2(bw10,ones(1,30),'same');
else
bw11=conv2(bw10,ones(1,75),'same');  
end
bw11=bwareaopen(bw11,200);%������
% imshow(bw11)

conn=bwconncomp(bw11);%�ҵ���ɫ��ͨ��
clist=conn.PixelIdxList;
L=0;
for i=1:length(clist)
    [m,n]=ind2sub(size(bw11),clist{i});
    mT=max(n)+(size(bw11,1)-min(m)); %����ӽ����Ͻǵ���ͨ��
    if mT>L
        L=mT;
        cidx=i;
        f_m=m;
        f_n=n;
    end
end
im7=im6(min(f_m):max(f_m),min(f_n):max(f_n),:);

ocrResults=ocr(im6);%ocrʶ��im6�����ַ�
wpos=ocrResults.WordBoundingBoxes;%�ҵ�ʶ���λ��
[~,idx]=min(abs(min(f_m)-wpos(:,2))+abs(min(f_n)-wpos(:,1)));%�ҵ�ocr����ӽ���λ��
figure
imshow(im7)
try
title(ocrResults.Words{idx})
fileID = fopen('result.txt','w');
fprintf(fileID,['ʶ��Ľ��Ϊ�� ',ocrResults.Words{idx}]);
fclose(fileID);
catch
end
else
    title('û����Ч��ǩ')
end
% title(ocrResults.Words)
%% ��ά��
qr=imread('qrcode.jpg');
qr=rgb2gray(qr);
points = detectSURFFeatures(qr); %���surf������
[features, valid_points] = extractFeatures(qr, points); %��ȡ����������(n*64)������������
idx = matchFeatures(features,features2);
figure
imshow(im)
hold on
plot(valid_points2(idx(:,2)),'showOrientation',true);
title('��ά����')
location=valid_points2(idx(:,2)).Location; %�������λ��
if numel(location)>2
s1 = std2(location(:,1)); %����xy������
s2=std2(location(:,2)); 
%%
% if isSmall==1 %�����Сͼ ����С�˴���
% gausI=imgaussfilt(rgb2gray(im),5); %��˹����,�ҵ���ɫ����
% else
%     gausI=imgaussfilt(rgb2gray(im),15);
% end
bw=rgb2gray(im)<20; %��ֵɸѡ
if isSmall==1
bw2=bwmorph(bw,'dilate',15);
else
bw2=bwmorph(bw,'dilate',25); %���ʹ���
end
figure
imshow(bw2)


%%
conn=bwconncomp(bw2); %����ͨ����
clist=conn.PixelIdxList;
maxC=0;
for i=1:length(clist)
    if numel(clist{i})>1000
        [m,n]=ind2sub(size(bw),clist{i}); %�±�ת��Ϊ����
%         if std(m)<s2*2 && std(n)<s1*2
            %�ڵ�ǰ��ͨ�����Ƿ���surf������
            mpos=find(location(:,1)>min(n) & location(:,1)<max(n)); 
            npos=find(location(:,2)>min(m) & location(:,2)<max(m));
            c=intersect(mpos,npos);%�ҵ�ͬʱ����mn�����������
            if numel(c)>maxC%ͳ����������������
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
    title('δ��⵽��ά��')
end
toc