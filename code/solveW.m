function W = solveW(Z,Y3,F,T,beta,mu)

n=size(Z,1);
R=zeros(n,n);
for i=1:n
    for j=1:n
        R(i,j)=1/2*norm((F(i,:)-F(j,:)),2)^2;
    end
end

temp = Z+Y3/mu;
[m,n] = size(Z);
temp_W = zeros(m,n);
for i=1:m
    for j=1:n
        if temp(i,j)>beta/(2*mu)*R(i,j)
            temp_W(i,j) = temp(i,j) - beta/mu*R(i,j);
        elseif temp(i,j)<-beta/(2*mu)*R(i,j)
            temp_W(i,j) = temp(i,j) + beta/mu*R(i,j);
        end
    end
end
for i=1:m
    [~,index] = sort(abs(temp_W(:,i)));
     temp_W(index(1:(n-T)),i) = 0;
end
W = temp_W;
end