function acc = cal_acc(Y,gt,train_index,test_index)
Y = Y*diag(sparse(1./sqrt(sum(Y.^2))));
Z_gallery = Y(:,train_index);
Z_probe = Y(:,test_index);
Y_gallery = gt(train_index);
Y_probe = gt(test_index);
knn_model = fitcknn(Z_gallery',Y_gallery,'NumNeighbors',1);
Y_probe_pseudo = knn_model.predict(Z_probe');
acc = length(find(Y_probe_pseudo==Y_probe))/length(Y_probe);
end