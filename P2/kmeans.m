  function[Y, u]=kmeans(X, k, lamda)  

    [a,b] = size(X);    
    Y = zeros(1, b);     
    
    u = X(:,randperm(b,k));
     
     while 1
        dis_x_u = zeros(1,k);
        num = zeros(1,k);
        new_u = u;
        
        for m = 1:b
            for n = 1:k
                dis_x_u(n) = (X(1:2,m) - u(1:2,n))'* (X(1:2,m) - u(1:2,n)) + lamda * (X(3:4,m) - u(3:4,n))'* (X(3:4,m) - u(3:4,n));
            end
            [~, temp] = min(dis_x_u);
            Y(1,m) = temp;
        end
        K=0;
        for n = 1:k
            for m = 1:b
                if Y(1,m) == n
                    new_u(:,n) = new_u(:,n) + X(:,m);
                    num(n) = num(n) + 1;
                end
            end
            new_u(:,n) = new_u(:,n) / num(n);
            if norm(new_u(:,n) - u(:,n)) < 0.001
                K = K+1;
            end
        end
        if K == k
            break;
        else
            u = new_u;
        end
     end
                               
end

