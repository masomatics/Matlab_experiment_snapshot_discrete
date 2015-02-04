currentpart = 8

thetahistory = []
energyhistory = [] 
for k = 1:currentpart 
hoge =load(['theta_history_minus13_0_0_0_0_uneven_part', num2str(k)]) 
thetahistory = [thetahistory; hoge.thetahistory];
end

for k = 2:currentpart 
hoge =load(['energy_history_minus13_0_0_0_0_uneven_part', num2str(k)]) 
energyhistory = [energyhistory; hoge.energyhistory'];
end


[numiter, numpar] = size(thetahistory)


theta_a = [   0    0.5000   25.0000    6.0000    0.2000]
distancehistory = abs(thetahistory - repmat(theta_a, [numiter, 1])).^2;


figure(99); 
set(gcf,'Position',get(0,'ScreenSize'))

for p = 1:numpar
subplot(2, 3, p)
plot(thetahistory(:,p))
        title(['theta' num2str(p)]);         
end 

subplot(2,3,6)
plot(energyhistory);
title('Energy history')


figure(98); 
set(gcf,'Position',get(0,'ScreenSize'))

for p = 1:numpar
subplot(2, 3, p)
plot(distancehistory(:,p))
        title(['theta' num2str(p) 'distance history']);         
end 
subplot(2,3,6)
plot(sum(distancehistory,2));
title('L2 Distance history')