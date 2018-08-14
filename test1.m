[m,n]=find(bw);
bw=bw4;
exit=0;
for i=1:size(bw,1) %’“◊Û…œΩ«∂•µ„
    for j=1:i
        try
            if bw(j,i+1-j)==1
                ul=[i-j,j];
                exit=1;
            end
        catch
        end
        if exit==1
            break
        end
    end
    if exit==1
        break
    end
end
%
temp_bw=fliplr(bw);
exit=0;
imshow(temp_bw);hold on
for i=1:size(bw,2) %’“”“…œΩ«∂•µ„
    for j=1:i
    if temp_bw(j,i+1-j)==1
        scatter(j,i+1-j)
        ur=[size(bw,2)-i+j,j];
        exit=1;
        break
    end
    end
    if exit==1
        break
    end
end
%
temp_bw=flipud(bw);
exit=0;
for i=1:size(bw,1) %’“◊Ûœ¬Ω«∂•µ„
    for j=1:i
    if temp_bw(j,i+1-j)==1
        ll=[i-j,size(bw,1)-j];
        exit=1;
        break
    end
    end
    if exit==1
        break
    end
end
%
% temp_bw=flipud(bw);
temp_bw=rot90(bw,2);
exit=0;
for i=1:size(bw,1) %’“”“œ¬Ω«∂•µ„
    for j=1:i
    if temp_bw(j,i+1-j)==1
        lr=[size(bw,2)-i+j,size(bw,1)-j];
        exit=1;
        break
    end
    end
    if exit==1
        break
    end
end
%%
imshow(bw)
hold on 
scatter(ur(1),ur(2))
scatter(ul(1),ul(2))
scatter(ll(1),ll(2))
scatter(lr(1),lr(2))
%%
new_ur=[ul(1)+pdist2(ur,ul),ul(2)];
% new_ll=[ul(1),ul(2)+pdist2(ul,ll)];
new_lr=[new_ur(1),new_ur(2)+pdist2(new_ur,lr)];
tf=([[ul',new_ur',new_lr'];[1 1 1]])*inv([[ul',ur',lr'];[1 1 1]]);
%%
new_ur=[lr(1),ul(2)];
new_ll=[ul(1),lr(2)];
tf=([[new_ur;new_ll;lr]';[1 1 1]])/([[ur;ll;lr]';[1 1 1]])
%%
% tf=([[ul',new_ur',new_lr'];[1 1 1]])*inv([[ul',ur',lr'];[1 1 1]]);
% tf=([[new_ur',new_ll',lr'];[1 1 1]])*inv([[ur',ll',lr'];[1 1 1]]);
tf=tf';
tf(7:8)=0;
tf(9)=1
%%
% tform = affine2d([1 0 0; .5 1 0; 1 0 1]);
tform = affine2d(tf);
J = imwarp(im2,tform);
figure
subplot(121)
imshow(im2)
subplot(122)
imshow(J)