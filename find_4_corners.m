function [ul,ur,lr,ll]=find_4_corners(bw)
b2=imrotate(bw,45,'loose');
% imshow(b2)
[m,n]=find(b2);
min_n=find(n==(min(n)));
min_m=find(m==(min(m)));
max_n=find(n==(max(n)));
max_m=find(m==(max(m)));
b2(m(min_n(1))-1:m(min_n(1))+1,n(min_n(1))-1:n(min_n(1))+1)=2; %左上，
b2(m(min_m(1))-1:m(min_m(1))+1,n(min_m(1))-1:n(min_m(1))+1)=3; %右上
b2(m(max_n(1))-1:m(max_n(1))+1,n(max_n(1))-1:n(max_n(1))+1)=4; %右下
b2(m(max_m(1))-1:m(max_m(1))+1,n(max_m(1))-1:n(max_m(1))+1)=5; %左下
% imshow(b2)
% hold on
% scatter(n(min_n(1)),m(min_n(1)))
%%
b3=imrotate(b2,-45,'loose');
m2=round((size(b3,1)-size(bw,1))/2);
n2=round((size(b3,2)-size(bw,2))/2);
b3=b3(m2:m2+size(bw,1)-1,n2:n2+size(bw,2)-1);
[ul(2),ul(1)]=find(b3==2,1); %左上
[ur(2),ur(1)]=find(b3==3,1); %右上
[lr(2),lr(1)]=find(b3==4,1); %右下
[ll(2),ll(1)]=find(b3==5,1); %左下
%%
% imshow(bw)
% hold on 
% scatter(ur(1),ur(2))
% scatter(ul(1),ul(2))
% scatter(ll(1),ll(2))
% scatter(lr(1),lr(2))