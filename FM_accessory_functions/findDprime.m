function d_prime=findDprime(responseMat)
   
    % The function findDprime (thanks to Dr. Ido Maor) receives a 3D matrix
    % responseMAt with dimensions - N x S x rep, where:
    % N - num of units in a regions
    % S - num of stimuli
    % rep - number of trials/repetitions of each stimulus
    
    % The function uses this matrix to calculate pairwise d primes for the
    % population between each 2 stimuli. 
    % The output "d_prime" is a SxS matrix symmetric matrix showing the
    % pairwise d primes between stimuli
    
    cell_num=size(responseMat,1);
    stim_num=size(responseMat,2);
    rep_num=size(responseMat,3);

    for p=1:stim_num
            
        for q=1:stim_num

            avg_resp_p=squeeze(mean(responseMat(:,p,:),3));
            avg_resp_q=squeeze(mean(responseMat(:,q,:),3));
            distance_pq=sqrt(sum((avg_resp_p-avg_resp_q).^2));

            for i=1:rep_num

                resp_p(i)=sqrt(sum((squeeze(responseMat(:,p,i))-avg_resp_p).^2));
                resp_q(i)=sqrt(sum((squeeze(responseMat(:,q,i))-avg_resp_q).^2));

            end
            innerdist_p=mean(resp_p);
            innerdist_q=mean(resp_q);
            d_prime(p,q)=distance_pq/mean([innerdist_p innerdist_q]);
            
        end
        
    end
    
end