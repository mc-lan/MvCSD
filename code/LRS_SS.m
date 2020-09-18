function [ Z,F,E] = LRS_SS(Xl,Xu,Yl,Yu,P,F,alpha,lambda,k,T,DEBUG)
% The code is written by Mengcheng Lan,
% if you have any problems, please don't hesitate to contact me: lanmengchengds@gmail.com 
% If you find the code is useful, please cite the following reference:
% Min Meng, Mengcheng Lan, Jun Yu, Jigang Wu, 
% Multi-View Consensus Structure Discovery [J], 
% IEEE Transactions on Cybernetics, 2020.

maxIter = 60;
rho = 1.2;
mu = 1e-2;
max_mu = 1e6;
tol = 1e-6;
nl = size(Yl,1);
nu = size(Yu,1);
numOfView = length(Xl);
numOfSamples = nl+nu;
%% Initializing optimization variables
% intializing
for i=1:numOfView
    X{i} = [Xl{i} Xu{i}];
    X{i} = X{i}./(repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1));
    E{i} = zeros(k, numOfSamples);
    Y1{i} = zeros(k, numOfSamples);
end
Z = zeros(numOfSamples, numOfSamples);
W = zeros(numOfSamples, numOfSamples);
J = zeros(numOfSamples, numOfSamples);
Y2 = zeros(numOfSamples, numOfSamples);
Y3 = zeros(numOfSamples, numOfSamples);

%% Start main loop
iter = 0;
if DEBUG
    disp(['initial,rank(Z)=' num2str(rank(Z))]);
end
while iter<maxIter
    iter = iter + 1;
    for i=1:numOfView
        PX{i} = P{i}'*X{i};
    end
    
    %solving J
    M1 = Z+Y2/mu;
    [U1, S1, V1] = svd((M1+eps),'econ');
    S1 = diag(S1);
    svp_J = length(find(S1>1/mu));
    if svp_J>=1
        S1 = S1(1:svp_J)-1/mu;
    else
        svp_J = 1;
        S1 = 0;
    end
    J = U1(:, 1:svp_J)*diag(S1)*V1(:, 1:svp_J)';
    
    %solving Z
    M3 = J+W-(Y2+Y3)/mu;
    W1 = zeros(numOfSamples,numOfSamples,numOfView);
    W2 = zeros(numOfSamples,numOfSamples,numOfView);
    for i=1:numOfView
        W1(:,:,i) = PX{i}'*(PX{i}-E{i}+Y1{i}/mu);
        W2(:,:,i) = PX{i}'*PX{i};
    end
    WW1 = sum(W1,3)+M3;
    WW2 = sum(W2,3)+2*eye(size(Z,2));
    Z = WW2\WW1;
    
    
    %solving W
    W = solveW(Z,Y3,F,T,lambda,mu);
    %solving F
    F = solveF(W,Yl);
    
    for i=1:numOfView
        %solving E
        temp_E = PX{i}-PX{i}*Z+Y1{i}/mu;
        for j = 1 : size(temp_E, 2)
            if norm(temp_E(:,j), 2) > alpha/mu
                E{i}(:,j) =( norm(temp_E(:,j), 2) - alpha/mu)/norm(temp_E(:,j),2)*temp_E(:,j);
            else
                E{i}(:,j) = 0;
            end
        end
    end
    
    dY1 = zeros(k,numOfSamples,numOfView);
    recErr1 = zeros(1,numOfView);
    for i=1:numOfView
        dY1(:,:,i) = PX{i}-PX{i}*Z-E{i};
        recErr1(i) = norm(dY1(:,:,i),'inf');
    end
    recerr1 = max(recErr1);
    dY2 =  Z - J;
    recErr2 = norm(dY2,'inf');
    dY3 =  Z - W;
    recErr3 = norm(dY3,'inf');
    recErr = max(recErr3,max(recerr1, recErr2));
    convergenced = recErr <tol||iter>maxIter;
    
    if DEBUG
        if iter==1 || mod(iter,5)==0 || convergenced
            disp(['iter ' num2str(iter) ',mu=' num2str(mu) ...
                ',rank(Z)=' num2str(rank(Z)) ...
                ',recErr=' num2str(recErr)]);
        end
    end
    if convergenced
        break;
    else
        for i=1:numOfView
            Y1{i} = Y1{i}+mu*dY1(:,:,i);
        end
        Y2 = Y2 + mu*dY2;
        Y3 = Y3 + mu*dY3;
        mu = min(max_mu, mu*rho);
    end
end
end
