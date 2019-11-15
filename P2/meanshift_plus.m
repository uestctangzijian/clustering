function [labels,clusters] = meanshift_plus(X, hc1, hc2, hp1, hp2, h)

[a,b] = size(X);
% I = eye(a);
numCluster = 0; 
clusters = []; 
threshold = 1e-3*h;
labels = zeros(1,b); 
sign = zeros(1,b); 

for i = 1:b
    if sign(i) == 1
        continue; 
    end
    t = 0;
    u0 = X(:,i);
%     sigma = (h^2) * I; 
    k = []; 
    x_up = zeros(a,1);
    x_down = zeros(a,1);
    while true
        if t == 0 
            old_miu = u0;
        else
            old_miu = new_miu;
        end
        
        d = sum((repmat(old_miu,1,b) - X).^2);
        k = find(d < (h^2)); 
        for ii = unique(k)
            dis_c1 = (X(1,ii)-old_miu(1))'* (X(1,ii)-old_miu(1));
            dis_c2 = (X(2,ii)-old_miu(2))'* (X(2,ii)-old_miu(2));
            dis_p1 = (X(3,ii)-old_miu(3))'* (X(3,ii)-old_miu(3));
            dis_p2 = (X(4,ii)-old_miu(4))'* (X(4,ii)-old_miu(4));
            kernel = (1/(4*(pi^2)*(hp1^2)*(hp2^2)*(hc1^2)*(hc2^2)))*exp(-0.5*(1/hc1^2)*dis_c1-0.5*(1/hc2^2)*dis_c2-0.5*(1/hp1^2)*dis_p1-0.5*(1/hp2^2)*dis_p2);
            x_up = x_up + X(:,ii) * kernel;
            x_down = x_down + kernel;
        end
        new_miu = x_up./x_down;
        
         if norm(new_miu - old_miu) < threshold
            
            merge = 0; 
            for iii = 1:numCluster
                dist2other = norm(new_miu-clusters(:,iii));
                if dist2other < h/2  
                    merge = iii;
                    break
                end  
            end
                  
            if merge > 0    
                clusters(:,merge) = 0.5*(new_miu+clusters(:,merge)); 
                labels(i) = merge;    
                for s = k
                    if sign(s) == 0
                        sign(s) = 1;
                        labels(s) = merge;
                    end
                end
            else    
                numCluster = numCluster+1; 
                clusters(:,numCluster) = new_miu;     
                labels(i) = numCluster;
                for s = k
                    sign(s) = 1;
                    labels(s) = numCluster;
                end
            end
            break;
        end
        t = t + 1; 
    end
        
end

