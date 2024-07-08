clear 

R =1; 
Pvec = [10 20 30 40];
ct=1000000000;
for i = 1 : length(Pvec)  
    sum1 = 0; sum2=0; sum3=0;
    Pbudget = 10^(Pvec(i)/10);
    parfor ict = 1 : ct
        M = 2;
        
        hall =  complex(randn(M,M)/sqrt(2),randn(M,M)/sqrt(2));
        hall = abs(hall).^2;
         % hall = [1 4 0.1 2 0.5]';%[1 2 3 4 5 ]'; %all the channels
         % hall = kron(hall,ones(1,M))
    
        %OMA
        Poma=0; % max power needed by OMA
        for m =1 : M
            Poma(m) = (exp(R)-1)/hall(m,m);
        end
    
        % %hybrid NOMA 
        P = zeros(M,M); 
        P(1,1) = (exp(R)-1)/hall(1,1); %User 1's transmit power 
        for m = 2: M 

            %%% method 2: using fmincon as the alternative 
            nonlcon = @mycons;%(x,N,tm,m,P,hall);
            options = optimoptions('fmincon','Display', 'off','MaxFunctionEvaluations', 300000); %display off
            x0 = zeros(m,1);
            A = []; % No other constraints
            b = [];
            Aeq = [];
            beq = [];
            lb = [];
            ub = [];
            x=[];
            x = fmincon(@(x) sum(x'),x0,A,b,Aeq,beq,lb,ub,@(x) mycons(x,m,P,hall,R),options);
            P(m,1:m)=x';

        end 
        

        %%%
        %Pana = [0 0; exp(R)/sqrt(hall(2,1)*hall(2,2))-exp(R)/hall(2,1) exp(R)/sqrt(hall(2,1)*hall(2,2))-1/hall(2,2) ]
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
        end
        
        %%%%
        Pnoma = sum(x);
        Pomaa = Poma(2);
        Pnomaa = sum(y);
        if Pnoma>Pbudget
            sum1 = sum1 +1;
        end
        if Pomaa>Pbudget
            sum2 = sum2 +1;
        end
        if Pnomaa>Pbudget
            sum3 = sum3 +1;
        end
    end
    Pout_noma(i) = sum1/ct;
    Pout_oma(i) = sum2/ct;
    Pout_nomaa(i) = sum3/ct;
    Pout_omaa(i) = (exp(R)-1)/Pbudget;
end
 
semilogy(Pvec,Pout_oma, Pvec, Pout_omaa,  Pvec, Pout_noma, Pvec, Pout_nomaa,...
    Pvec, 1./10.^(Pvec/10),Pvec, 1./10.^(2*Pvec/10))
%semilogy(  Pvec, Pout_omaa,  Pvec, Pout_nomaa,...
%    Pvec, 1./10.^(Pvec/10),Pvec, 10./10.^(2*Pvec/10))

function [c,ceq] = mycons(x,m,P,hall,R)

c(1) = R;
for i = 1: m  
    hmn = hall(m,i);
    c(1) = c(1) - log(1 + hmn*x(i)/(1+sum(hall(i:m-1,i).*P(i:m-1,i)))) ;
    c(i+1) = -x(i);
end
    ceq = [];
 
end

 