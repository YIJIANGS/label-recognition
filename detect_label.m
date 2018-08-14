clear;clc;close all;
im=imread('i2.jpg'); %��ȡͼƬ
figure
imshow(im)
title('ԭͼ')
%�ֱ���ȡrgb��ɫͨ��
r=im(:,:,1); 
g=im(:,:,2);
b=im(:,:,3);
%% �ҳ���ɫ��ǩ��Χ���ֶ����ԣ�
r1=r<50;
g1=g>50 & g<160;
b1=b>70;
bw=r1.*g1.*b1; %ͬʱ������ɫ��Χ�ĵ�
%% ��ȡ����ʶ
bw2=conv2(bw,ones(2),'same'); %��������ʶ��ʶ��Χ�������൱������
bw3=imfill(bw2,'holes'); %�������ͼ
conn=bwconncomp(bw3); %�ҵ���ͨ����
plist=conn.PixelIdxList; %������ͨ�����±�
pnum=cellfun(@length,plist); %������ͨ��������ͨ������
[~,idx]=max(pnum); %�ҳ������ͨ��������
%���¹�����ɫͨ�� ���зǱ�ʶ����Ϳ��
r3=zeros(size(r1)); g3=zeros(size(r1)); b3=zeros(size(r1));
% r3=nan(size(r1)); g3=nan(size(r1)); b3=nan(size(r1));
%���е�ɫͨ��������ɫ
r3(plist{idx})=r(plist{idx});
g3(plist{idx})=g(plist{idx});
b3(plist{idx})=b(plist{idx});
im2=uint8(cat(3,r3,g3,b3)); %��ɫͨ���ں�
[m,n]=find(rgb2gray(im2)); %�ҵ���ɫ��������
im2=im(min(m):max(m),min(n):max(n),:); %��ȡ��ɫ����ͼ
figure
imshow(im2) %��ʾ����ʶͼ��
title('���ʶͼ')
%% ���ʶ��ת����
bw4=bw3(min(m):max(m),min(n):max(n),:); %��ԭ�����ֵͼ�н�ȡ���ʶ����
% bw4=zeros(size(bw3));
% bw4(plist{idx})=1;
bw4 = bwconvhull(bw4); %ͨ����Ե����϶����
minL=2000;
%����һ���Ƕȷ�Χ ��ÿ0.5�ĽǶȽ���ͼ����ת��ÿ����ת��ͳ����ˮƽ�ʹ�ֱ��������ص�
%������������ռ�����ٵĽǶ�
for i=-30:.5:30
    tempI=imrotate(bw4,i,'loose');
    [m,n]=find(tempI);
    L=length(unique(n))+length(unique(m));
    if L<minL
        minL=L;
        ang=i;
    end
end
im3=imrotate(im2,ang,'loose'); %ȷ���ǶȺ��ͼ�������ת
% imshow(im3)
%% ��ȡС���ʶ
%�ֱ���ȡrgb��ɫͨ��
r=im3(:,:,1); 
g=im3(:,:,2);
b=im3(:,:,3);
%�ҳ���ɫ��ǩ��Χ���ֶ����ԣ�
r1=r>185;
g1=g>180;
b1=b>155;
bw5=r1.*g1.*b1; %ͬʱ������ɫ��Χ�ĵ�
bw6=conv2(bw5,ones(2),'same'); %��������ʶ��ʶ��Χ�������൱������
bw7=imfill(bw6,'holes'); %����϶
conn=bwconncomp(bw7); %�ҵ���ͨ����
plist=conn.PixelIdxList; %������ͨ�����±�
pnum=cellfun(@length,plist); %������ͨ��������ͨ������
[~,idx]=sort(pnum,'descend'); %�ҳ����2����ͨ��������
%%
%�ֱ��2��������ͨ������±�ת��Ϊ����
[x,y]=ind2sub(size(bw7),plist{idx(1)}); 
[x1,y1]=ind2sub(size(bw7),plist{idx(2)});
figure
imshow(im3)
hold on
%�ֱ��ҵ�2������ı߽�
k = boundary(x,y,.02); %��ϱ߽�
plot(y(k),x(k),'r') 
k=boundary(x1,y1,.02);
plot(y1(k),x1(k),'r')
if mean(x)<mean(x1) %ͨ��yֵƽ��ֵ�жϱ�ǩ��λ��
    idx1=idx(1);
    idx2=idx(2);
else
    idx1=idx(2);
    idx2=idx(1);
end
title('���ʶ��ת����С��ʶ��λ')
%% С��ʶ����
bw8=zeros(size(bw7)); %����һ��ȫ�ھ���
bw8(plist{idx1})=1;
bw9 = bwconvhull(bw8); %�߽����
[m,n]=find(bw9);
%���¹�����ɫͨ�� ���зǱ�ʶ����Ϳ��
r3=zeros(size(r1)); g3=zeros(size(r1)); b3=zeros(size(r1));
%���е�ɫͨ��������ɫ
r3(bw9)=r(bw9);
g3(bw9)=g(bw9);
b3(bw9)=b(bw9);
im4=uint8(cat(3,r3,g3,b3)); %��ɫͨ���ں� ��תΪuint8��
im4=im3(min(m):max(m),min(n):max(n),:);
% imshow(im4) %��ʾС���ʶͼ��
%% С��ʶ��ת����
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

%% ��Чͼ����
bw8=zeros(size(bw7));
bw8(plist{idx2})=1;
bw9 = bwconvhull(bw8);
[m,n]=find(bw9);
%���¹�����ɫͨ�� ���зǱ�ʶ����Ϳ��
r3=zeros(size(r)); g3=zeros(size(r)); b3=zeros(size(r));
%���е�ɫͨ��������ɫ
r3(bw9)=r(bw9);
g3(bw9)=g(bw9);
b3(bw9)=b(bw9);
im4=uint8(cat(3,r3,g3,b3)); %��ɫͨ���ں� ��תΪuint8��
im4=im3(min(m):max(m),min(n):max(n),:);
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
im6=imrotate(im4,ang,'crop');
figure
imshow(im6)
title('��Чͼ')
%% ʶ����Чͼ
gray=rgb2gray(im6);
bw10=gray<40 & gray>0;
bw11=conv2(bw10,ones(1,15),'same');
bw11=bwareaopen(bw11,200);
%%
conn=bwconncomp(bw11);
clist=conn.PixelIdxList;
L=0;
for i=1:length(clist)
    [m,n]=ind2sub(size(bw11),clist{i});
    mT=max(n)^2+(size(bw11,1)-min(m))^2; %����ӽ����Ͻǵ���ͨ��
    if mT>L
        L=mT;
        cidx=i;
        f_m=m;
        f_n=n;
    end
end
try
im7=im6(min(f_m)-10:max(f_m)+10,min(f_n)-10:max(f_n)+10,:);
catch
    im7=im6(min(f_m):max(f_m),min(f_n):max(f_n),:);
end
figure
imshow(im7)
ocrResults = ocr(im7);
title(ocrResults.Words)