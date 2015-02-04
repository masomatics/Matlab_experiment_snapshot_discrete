alpha = [0.005, 0.2, 0.3];
init = [0; 100]; 
tend = 1000;
rnsource1 = rand(1,50000);
rnsource2 = rand(1,50000);


[t, x] = SI_model_Arai(init, tend, alpha, rnsource1, rnsource2);
plot(t,x)
