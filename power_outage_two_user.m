clear all

R =2; 
Pvec = [0 :5: 30];
M=2;
ct=5000000;
parfor i = 1 : length(Pvec)  
    sum1 = 0; sum2=0; sum3=0; sumh=0; sum4=0; sum5=0;
    Pbudget = 10^(Pvec(i)/10);
    for ict = 1 : ct     
        hall =  complex(randn(M,M)/sqrt(2),randn(M,M)/sqrt(2));
        hall = abs(hall).^2; 
    
        %OMA
        Poma=0; % max power needed by OMA
        for m =1 : M
            Poma(m) = (exp(R)-1)/hall(m,m);
        end        

        %%%
        crit = hall(2,2)/hall(2,1)-exp(-2*R);
        y=zeros(2,1);
        if crit >=1
            %fprintf('Pure OMA')
            lambda1 = exp(R)/hall(2,2);
            lambda2 = 1 - hall(2,1)/hall(2,2);
            lambda3 = 0;
            y(1) = 0;
            y(2) = (exp(R)-1)/hall(2,2);
        elseif crit<= exp(-2*R)
            %fprintf('Pure NOMA')
            lambda1 = exp(2*R)/hall(2,1);
            lambda2 = 0;
            lambda3 = 1 - exp(2*R)*hall(2,2)/hall(2,1);
            y(1) = exp(R)*(exp(R)-1)/hall(2,1);
            y(2) = 0;
        else
            %fprintf('Hybrid NOMA')
            lambda2 = 0;
            lambda3 = 0;
            lambda1 = exp(R)/sqrt(hall(2,1)*hall(2,2));
            y(1) = lambda1 - exp(R)/hall(2,1);
            y(2) = lambda1 - 1/hall(2,2);
            c=0;
        end
        
        %analysis
        if Poma(2)>Pbudget
            sum1 = sum1 +1;
        end
        if sum(y)>Pbudget
            sum2 = sum2 +1;
        end

         
    end 
     Ppoma(i) = sum1/ct;
     Psim_all(i) = sum2/ct; 

    %analysis
    stepy = exp(R)/1000;
    sumy = 0; sumy2=0;
    yvec = [1:stepy: exp(R)];
    for iy = 1 : length(yvec)
        y = yvec(iy);
        g2y = (2*exp(2*R)+2*exp(R)*sqrt(exp(2*R)-exp(R)*y)-exp(R)*y)/y^2;
        g1y = (2*exp(2*R)-2*exp(R)*sqrt(exp(2*R)-exp(R)*y)-exp(R)*y)/y^2;
        sumy = sumy + (min(exp(2*R),g2y)-max(1,g1y))*(y-1)*stepy;
    end
    hybrid = sumy/Pbudget^2;
    pureoma = (exp(R)-1)^2/Pbudget^2/2; 
    Pana(i) = 2*pureoma + hybrid;

    Pomaa(i) = 1-exp(-(exp(R)-1)/Pbudget);
end
 
semilogy(Pvec, Ppoma,  Pvec, Psim_all , Pvec, Pomaa,  Pvec, Pana) 
 

 