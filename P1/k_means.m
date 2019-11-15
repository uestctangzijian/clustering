  function[Y, u]=k_means(X, k, m)  

    [a,b] = size(X);    
    Y = zeros(1, b);     
    u = zeros(a, k);
    
    if m == 1
        u(:,1)=[10,-10];
        u(:,2)=[-10,-10];
        u(:,3)=[-10,10];
        u(:,4)=[10,10];
    end
    
    if m == 2
        u(:,1)=[-5,-0];
        u(:,2)=[-1,-0];
        u(:,3)=[1,0];
        u(:,4)=[5,0];
    end
    
    if m == 3
        u(:,1)=[-10,0];
        u(:,2)=[-1,0];
        u(:,3)=[1,0];
        u(:,4)=[10,0];
    end   
     
     while 1
        dis_x_u = zeros(1,k);
        num = zeros(1,k);
        new_u = u;
        
        for m = 1:b
            for n = 1:k
                dis_x_u(n) = dot((X(:,m) - u(:,n)), (X(:,m) - u(:,n)));
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
            if norm(new_u(:,n) - u(:,n)) < 0.0001
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
