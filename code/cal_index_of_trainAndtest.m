function [train_index,test_index] = cal_index_of_trainAndtest(Y,per)
classes = length(unique(Y));
index_all = 1:length(Y);

toll = 0;
train_index = [];
for i=1:classes
    num0fclass(i) = length(find(Y==i));
    if i==1
        toll = 0;
    else
        toll = toll + num0fclass(i-1);
    end
    index1 = randperm(num0fclass(i));
    numOflabeled = ceil(num0fclass(i)*per);
    train_index_i = toll + index1(1:numOflabeled);
    train_index = [train_index train_index_i];
end
for i=1:length(train_index)
    index_all(index_all==train_index(i)) = [];
end
test_index = index_all;
end