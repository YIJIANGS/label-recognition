%% im 原图 bw 二值图 ul,ur,lr,ll 4个顶点（左上，右上，右下，左下)
function final_im=perspective_transform(im,bw,ul,ur,lr,ll)
bw=bwconvhull(bw);
% img=bw;
dot=[ul;ur;lr;ll];
new_ur=[ur(1),ul(2)];
new_ll=[ul(1),ll(2)];
new_lr=[ur(1),ll(2)];
dot2=[ul;new_ur;new_lr;new_ll];

y=[dot(1,1) dot(2,1) dot(3,1) dot(4,1)];        %四个原顶点
x=[dot(1,2) dot(2,2) dot(3,2) dot(4,2)];

%这里是新的顶点，我取的矩形,也可以做成其他的形状
%大可以原图像是矩形，新图像是从dot中取得的点组成的任意四边形.:)
Y=[dot2(1,1) dot2(2,1) dot2(3,1) dot2(4,1)];        %四个新顶点
X=[dot2(1,2) dot2(2,2) dot2(3,2) dot2(4,2)];

B=[X(1) Y(1) X(2) Y(2) X(3) Y(3) X(4) Y(4)]';   %变换后的四个顶点，方程右边的值
%联立解方程组，方程的系数
A=[x(1) y(1) 1 0 0 0 -X(1)*x(1) -X(1)*y(1);             
  0 0 0 x(1) y(1) 1 -Y(1)*x(1) -Y(1)*y(1);
   x(2) y(2) 1 0 0 0 -X(2)*x(2) -X(2)*y(2);
  0 0 0 x(2) y(2) 1 -Y(2)*x(2) -Y(2)*y(2);
   x(3) y(3) 1 0 0 0 -X(3)*x(3) -X(3)*y(3);
  0 0 0 x(3) y(3) 1 -Y(3)*x(3) -Y(3)*y(3);
   x(4) y(4) 1 0 0 0 -X(4)*x(4) -X(4)*y(4);
  0 0 0 x(4) y(4) 1 -Y(4)*x(4) -Y(4)*y(4)];

fa=inv(A)*B;        %用四点求得的方程的解，也是全局变换系数
a=fa(1);b=fa(2);c=fa(3);
d=fa(4);e=fa(5);f=fa(6);
m=fa(7);l=fa(8);

% new_x=(a*x+b*y+c)/(m*x+l*y+1);
% new_y=(d*x+e*y+f)/(m*x+l*y+1);

[M, N,~] = size(bw);
newI=zeros(M,N,3);
for i=1:M
    for j=1:N
        if bw(i,j)==1
            x=i;y=j;
            new_x=round((a*x+b*y+c)/(m*x+l*y+1));
            new_y=round((d*x+e*y+f)/(m*x+l*y+1));
            newI(new_x,new_y,:)=im(x,y,:);
        end
    end
end
%%
b2=newI(:,:,1)>0;
b3=bwconvhull(b2);
b4=b3-b2;
%%
[m,n]=find(b4);
for i=1:length(m)
    temp_m=newI(m(i)-1:m(i)+1,n(i)-1:n(i)+1,:);
    temp_m(temp_m==0)=nan;
    newI(m(i),n(i),:)=nanmean(nanmean(temp_m,1));
end
[m,n]=find(b3);
final_im=newI(min(m):max(m),min(n):max(n),:);
final_im=uint8(final_im);
