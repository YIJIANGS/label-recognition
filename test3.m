b2=imrotate(bw4,45,'loose');
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
%%
[m,n]=find(b2==4);
imshow(b2)
hold on
scatter(n,m)
%%
b3=imrotate(b2,-45,'loose');
m2=round((size(b3,1)-size(bw,1))/2);
n2=round((size(b3,2)-size(bw,2))/2);
b3=b3(m2:m2+size(bw,1)-1,n2:n2+size(bw,2)-1);
[m,n]=find(b3==4,1);
imshow(b3)
hold on
scatter(n,m)
%%
imshow(b3)
hold on 
scatter(ur(1),ur(2))
scatter(ul(1),ul(2))
scatter(ll(1),ll(2))
scatter(lr(1),lr(2))