clear 
close all
%load data_exhaustive.mat
R =1; 
Mvec = [2 : 10];
load("data_2user_search.mat")
     M = 2;
     %hall =  complex(randn(M,M)/sqrt(2),randn(M,M)/sqrt(2));
     %hall = abs(hall).^2;
     % hall = [1 4 0.1 2 0.5]';%[1 2 3 4 5 ]'; %all the channels
     % hall = kron(hall,ones(1,M))

    %OMA
    Poma=0; % max power needed by OMA
    for m =1 : M
        Poma(m) = (exp(R)-1)/hall(m,m);
    end

    %hybrid NOMA 
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

 

 x


stepsize= 0.01;%sum(Poma)/1000;
%exhaustive search
pvec = [0:stepsize:1];%sum(Poma)+2*stepsize];

for i1 = 1   
    
    for i2 = 1:length(pvec)
        for i3 = 1:length(pvec) 
                        p11 = Poma(1);p21 = pvec(i2);p22 = pvec(i3); 
                        R2 = log(1+p22*hall(2,2))+log(1+p21*hall(2,1)/(1+p11*hall(1,1))); 
                        if R2>=R %feasible
                            powersearch(i2,i3) = p21+p22;
                        else
                            powersearch(i2,i3) = inf;
                        end
                         
                  
        end
    end
end

mesh(pvec, pvec, powersearch)

hold

mesh(pvec, pvec, sum(x)*ones(length(pvec),length(pvec)))

function [c,ceq] = mycons(x,m,P,hall,R)

c(1) = R;
for i = 1: m  
    hmn = hall(m,i);
    c(1) = c(1) - log(1 + hmn*x(i)/(1+sum(hall(i:m-1,i).*P(i:m-1,i)))) ;
    c(i+1) = -x(i);
end
    ceq = [];
 
end

 