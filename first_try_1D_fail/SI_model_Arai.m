
function[time, state] = SI_model_Arai(init, tend, alpha, rnsource1, rnsource2)


a(1) = alpha(1)*init(1)*init(2);
a(2) = alpha(2)*init(1);
a(3) = alpha(3)*init(2);

dimen = length(init); 
statenow = init; 
tnow = 0; 
maxiter = 10000;
state = zeros(dimen, maxiter);
state(:,1) = init;  
time = zeros(1, maxiter);


cum = [a(1), a(1)+a(2)]/(sum(a));  

zeta = [ -1 -1 1; 1 1 -1];  

    for(k = 2: maxiter)    
        deltat = -log(rnsource1(k))/(sum(a));  
        cumnow = 1; 
        choice = 1;
        for(j = 1 : 2)  
            if(rnsource2(k) > cum(cumnow))
                choice = choice +1; 
            end
            cumnow = cumnow  + 1;
        end 
        tnow = tnow  + deltat;  
        time(k) = tnow;

        statenow = statenow + zeta(:,choice);
        state(:,k) = statenow; 

        a(1)= alpha(1)*statenow(1)*statenow(2);
        a(2) = alpha(2)*statenow(1);
        a(3) = alpha(3)*statenow(2);
        cum = [a(1), a(1)+a(2)]/(sum(a));  

        if(tnow > tend)
            break; 
        end 
             
    end 
  time  = time(1:k);   
  state = state(:, 1:k);
    




end 