function [Xl,Xu,Yl,Yu] = DataCreate(X,gt,NumOfLabled)
view = length(X);
%NumOfLabled = 2;
class = length(unique(gt));
for i =1:view
    Xl{i} = [];
    Xu{i} = [];
end
Yl = [];
Yu = [];
for c =1:class
    index1 = find(gt==c);
    nc = length(index1);
    index2 = randperm(nc);
    for i=1:view
        Xl{i} = [Xl{i} X{i}(:,index1(index2(1:NumOfLabled)))];
        Xu{i} = [Xu{i} X{i}(:,index1(index2(NumOfLabled+1:end)))];
    end
    Yl = [Yl;gt(index1(index2(1:NumOfLabled)))];
    Yu = [Yu;gt(index1(index2(NumOfLabled+1:end)))];
end
end