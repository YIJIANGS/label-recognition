function varargout = label_gui(varargin)
% LABEL_GUI MATLAB code for label_gui.fig
%      LABEL_GUI, by itself, creates a new LABEL_GUI or raises the existing
%      singleton*.
%
%      H = LABEL_GUI returns the handle to a new LABEL_GUI or the handle to
%      the existing singleton*.
%
%      LABEL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LABEL_GUI.M with the given input arguments.
%
%      LABEL_GUI('Property','Value',...) creates a new LABEL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before label_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to label_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help label_gui

% Last Modified by GUIDE v2.5 04-May-2018 23:43:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @label_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @label_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before label_gui is made visible.
function label_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to label_gui (see VARARGIN)

% Choose default command line output for label_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
axes(handles.axes1);set(gca,'XTick',[]);set(gca,'YTick',[]);
axes(handles.axes2);set(gca,'XTick',[]);set(gca,'YTick',[]);
axes(handles.axes5);set(gca,'XTick',[]);set(gca,'YTick',[]);
% UIWAIT makes label_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = label_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% ��ȡԭͼ
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im isSmall
axes(handles.axes1);cla(gca);title('')
axes(handles.axes2);cla(gca);title('')
axes(handles.axes5);cla(gca);title('')
set(handles.text5,'String','')
set(handles.text6,'String','')
set(handles.text7,'String','')

[file,path,~] = uigetfile('*');
im=imread([path,file]); %��ȡͼƬ
if numel(im)<1e7 %�ж�ͼ���Ƿ���Сͼ��ɸѡ��ͬ��ֵ��
    isSmall=1;
else
    isSmall=0;
end
axes(handles.axes1)
imshow(im)



% ��λ��Чͼ
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im features2 valid_points2 idx1 idx2 bw7 plist im3

load('feature.mat')
points2 = detectSURFFeatures(rgb2gray(im)); %���surf������
[features2, valid_points2] = extractFeatures(rgb2gray(im), points2); %��ȡ����������(n*64)������������
idx = matchFeatures(features,features2);
location=valid_points2(idx(:,2)).Location; %�������λ��
if numel(location)>4 %ȷ���б�ǩ
%�ֱ���ȡrgb��ɫͨ��
r=im(:,:,1); 
g=im(:,:,2);
b=im(:,:,3);
%�ҳ���ɫ��ǩ��Χ���ֶ����ԣ�
r1=r<50;
g1=g>50 & g<160;
b1=b>70;
bw=r1.*g1.*b1; %ͬʱ������ɫ��Χ�ĵ�    
    
% ��ȡ����ʶ
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
% ���ʶ��ת����
bw4=bw3(min(m):max(m),min(n):max(n),:); %��ԭ�����ֵͼ�н�ȡ���ʶ����
bw4 = bwconvhull(bw4); %ͨ����Ե����϶����
minL=99999;
%����һ���Ƕȷ�Χ ��ÿ0.5�ĽǶȽ���ͼ����ת��ÿ����ת��ͳ����ˮƽ�ʹ�ֱ��������ص�
%������������ռ�����ٵĽǶ�
for i=-30:.5:30
    tempI=imrotate(bw4,i,'loose');
    [m,n]=find(tempI);
    L=max(n)-min(n)+max(m)-min(m);
    if L<minL
        minL=L;
        ang=i;
    end
end
im3=imrotate(im2,ang,'loose'); %ȷ���ǶȺ��ͼ�������ת
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
bw6=conv2(bw5,ones(2),'same'); %��������ʶ��ʶ��Χ�������൱������
bw7=imfill(bw6,'holes'); %����϶
conn=bwconncomp(bw7); %�ҵ���ͨ����
plist=conn.PixelIdxList; %������ͨ�����±�
pnum=cellfun(@length,plist); %������ͨ��������ͨ������
[~,idx]=sort(pnum,'descend'); %�ҳ����2����ͨ��������

%�ֱ��2��������ͨ������±�ת��Ϊ����
[x,y]=ind2sub(size(bw7),plist{idx(1)}); 
[x1,y1]=ind2sub(size(bw7),plist{idx(2)});
%�ֱ��ҵ�2������ı߽�
k = boundary(x,y,.02); %�ҵ��߽��
plot(y(k),x(k),'r') %�±�ת����Ĺ������е��ã�����yx
k=boundary(x1,y1,.02);
plot(y1(k),x1(k),'r')
if mean(x)<mean(x1) %ͨ��yֵƽ��ֵ�жϱ�ǩ��λ��
    idx1=idx(1);%idx����������px���λ��
    idx2=idx(2);
else
    idx1=idx(2);
    idx2=idx(1);
end

axes(handles.axes1)
imshow(im3)


else
    axes(handles.axes1)
    title('û����Ч��ǩ')
end


% ��λ��ά��
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im features2 isSmall valid_points2 qrcode
load('features.mat')
idx = matchFeatures(qrfeatures,features2);
location=valid_points2(idx(:,2)).Location; %�������λ��
if numel(location)>10
s1 = std2(location(:,1)); %����xy������
s2=std2(location(:,2)); 

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

conn=bwconncomp(bw2); %����ͨ����
clist=conn.PixelIdxList;
maxC=0;
for i=1:length(clist)
    if numel(clist{i})>1000
        [m,n]=ind2sub(size(bw),clist{i}); %�±�ת��Ϊ����
        if std(m)<s2*2 && std(n)<s1*2
            %�ڵ�ǰ��ͨ�����Ƿ���surf������
            mpos=find(location(:,1)>min(n) & location(:,1)<max(n)); 
            npos=find(location(:,2)>min(m) & location(:,2)<max(m));
            c=intersect(mpos,npos);%�ҵ�ͬʱ����mn�����������
            if numel(c)>maxC%ͳ����������������
                maxC=numel(c);
                final_m=m;
                final_n=n;
            end
        end
    end
end

newI=im(min(final_m):max(final_m),min(final_n):max(final_n),:);

bw2=rgb2gray(newI)<40;
bw3=bwareaopen(bw2,150); %ȥ�봦��

minL=99999;
for i=-30:.5:30
    tempI=imrotate(bw3,i,'loose');
    [m,n]=find(tempI);
    L=max(n)-min(n)+max(m)-min(m);
    if L<minL
        minL=L;
        ang=i;
        tempBW=tempI;
    end
end
newI2=imrotate(newI,ang,'loose');
axes(handles.axes1)
imshow(newI2)
qrcode=1;
else
    qrcode=0;
    axes(handles.axes1)
    title('δ��⵽��ά��')
end



% ����ͼƬ1
function pushbutton4_Callback(hObject, eventdata, handles)
global idx1 bw7 plist im3
bw8=zeros(size(bw7)); %����һ��ȫ�ھ���
bw8(plist{idx1})=1;
bw9 = bwconvhull(bw8); %�߽����
[m,n]=find(bw9);
im4=im3(min(m):max(m),min(n):max(n),:);
minL=inf;
for i=-30:.5:30
    tempI=imrotate(bw9,i,'loose');
    [m,n]=find(tempI);
    L=max(n)-min(n)+max(m)-min(m);
    if L<minL
        minL=L;
        ang=i;
    end
end
im5=imrotate(im4,ang,'crop');
g=rgb2gray(im5);
g=imadjust(g);
b2=g>220;
b3=bwareaopen(b2,1000);
b3=bwconvhull(b3);
b4=imcomplement(b2).*b3;
axes(handles.axes2)
imshow(im5)
axis tight
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
        set(handles.text5,'String',['3-',num2str(find(orderId==Midx(1)))])
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
        set(handles.text5,'String',['5-',num2str(find(orderId==Midx(1)))])
end
title(label)

% ����ͼƬ2
function pushbutton6_Callback(hObject, eventdata, handles)
global bw7 im3 plist idx2 isSmall
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
im4=im3(min(m):max(m),min(n):max(n),:);
minL=99999;
for i=-30:.5:30
    tempI=imrotate(bw9,i,'loose');
    [m,n]=find(tempI);
    L=max(n)-min(n)+max(m)-min(m);
    if L<minL
        minL=L;
        ang=i;
    end
end
im6=imrotate(im4,ang,'crop');

% ʶ����Чͼ
gray=rgb2gray(im6);
bw10=gray<40 & gray>0;%��ֵ��
if isSmall==1
bw11=conv2(bw10,ones(1,30),'same');
else
bw11=conv2(bw10,ones(1,75),'same');  
end
bw11=bwareaopen(bw11,200);%������
conn=bwconncomp(bw11);%�ҵ���ɫ��ͨ��
clist=conn.PixelIdxList;
L=0;
for i=1:length(clist)
    [m,n]=ind2sub(size(bw11),clist{i});
    mT=max(n)+(size(bw11,1)-min(m)); %����ӽ����Ͻǵ���ͨ��
    if mT>L
        L=mT;
        f_m=m;
        f_n=n;
    end
end
im7=im6(min(f_m):max(f_m),min(f_n):max(f_n),:);

ocrResults=ocr(im6);%ocrʶ��im6�����ַ�
wpos=ocrResults.WordBoundingBoxes;%�ҵ�ʶ���λ��
[~,idx]=min(abs(min(f_m)-wpos(:,2))+abs(min(f_n)-wpos(:,1)));%�ҵ�ocr����ӽ���λ��
axes(handles.axes5)
imshow(im7)
title(ocrResults.Words{idx})
set(handles.text6,'String',ocrResults.Words{idx})


% ʶ��
function pushbutton7_Callback(hObject, eventdata, handles)
global qrcode
var1=get(handles.text5,'String');
var2=get(handles.text6,'String');
if strcmp(var1,'3-2') && strcmp(var2,'3.59') && qrcode==1
    set(handles.text7,'String','��ǩ��ȷ')
else
    set(handles.text7,'String','��ǩ����')
end
