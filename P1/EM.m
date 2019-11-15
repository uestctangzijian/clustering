function [miu, zi] = EM(X,k,init_miu,init_sigma,init_pi)

    [a,b] = size(X);

    zij = zeros(1,b);
    zi = zeros(1,b);

    miu = init_miu;
    pi = init_pi;
    sigma = init_sigma;

    while true
        old_miu = miu;
        old_sigma = sigma;
        old_pi = pi;
        
        %EMloop         
        for ii=1:b
            Pi_norm = zeros(1,k);
            sum_normK = 0;
            for kk = 1:k
                Pi_norm(:,kk) = old_pi(:,kk) * mvnpdf(X(:,ii),old_miu(:,kk),old_sigma(:,:,kk)); 
                sum_normK = sum_normK + Pi_norm(:,kk);
            end
            max_pn = max(Pi_norm);
            max_k = find(Pi_norm == max_pn);
            if max_pn == 0
                zij(1,ii) = 0;
            else
                zi(1,ii) = max_k;
                zij(1,ii) = max_pn./sum_normK;
            end
        end

        Nj = zeros(1,k);
        for kkk = 1:k
            j = find(zi == kkk); 
            for fj = j
                Nj(1,kkk) = Nj(1,kkk) + zij(1,fj); 
            end
        end

        PiJ = zeros(1,k);
        for kv = 1:k
            PiJ(1,kv) = Nj(1,kv) / b; 
        end

        zij_x = zeros(a,k); 
        for kvv = 1:k
            jj = find(zi == kvv);
            for fjj = jj
                zij_x(:,kvv) = zij_x(:,kvv) + zij(1,fjj) * X(:,fjj);
            end
            miu(:,kvv) = zij_x(:,kvv) / Nj(1,kvv);
        end

        zij_x_miu = zeros(a,a,k); 
        for kvi = 1:k
            jjj = find(zi == kvi);
            for fjjj = jjj
                zij_x_miu(:,:,kvi) = zij_x_miu(:,:,kvi) + zij(1,fjjj) * (X(:,fjjj)-miu(:,kvi)) * (X(:,fjjj)-miu(:,kvi))';
            end
            sigma(:,:,kvi) = zij_x_miu(:,:,kvi) / Nj(1,kvi);
        end
    
        
        if norm(miu-old_miu) < 0.00001
            break
        end
    end
    
end
        
        
