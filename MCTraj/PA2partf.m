function [] = PA2partf(nParticles)
% Baldeep Kooner 101004107
% ELEC 4700 MonteCarlo PA (2) Assignment
% nTime = number of time steps
% Ex = Electric Field in the x-direction
Ex = 5*10^-20; 
nTime = 100;
dt = 1; 
Ps = 0.05; %probability of scattering
%P = @(dx) 1 - exp(-dx/lambda); 
%P = 1 - 1.2*rand();
%x0 = 0; 
%v0 = 0;
t = 1:nTime;%time
dt = 1; 
x = zeros(nParticles, nTime); %distance
vx = zeros(nParticles, nTime); %velocity
for k = 1:nParticles
    
    for n = 2:nTime
      P = 1 - rand();
      if P <= 0.05
          vx(k, n) = 0;
      else
          vx(k, n) = vx(k, n-1) + (1.60217653e-19)*Ex*dt / (9.10938215e-31);
      end
      %x(n) = vx(n-1)*dt + ((1.60217653e-19)*Ex*dt^2) / (2*(9.10938215e-31));
      x(k, n) = x(k, n-1) + vx(k, n-1)*dt;
    end
end

vdx = zeros(nParticles, nTime);
for j = 1:nParticles
    vdx(j, :) = sum(vx(j, :)) / nTime;
end
vdxsum = 0;
for l = 1:nParticles
    vdxsum = vdxsum + sum(vdx(l, :));
end
totalVdx = vdxsum / nParticles;
%vdx = sum(vx)/nTime;
%plot(t, vx)
%xlabel('time')
%ylabel('velocity')
%hold on
%plot(t, x)

%comet(t, x)
for i = 1:nParticles
    comet(t, vx(i, :))
    xlabel('time')
    ylabel('velocity')
    hold on
end
title(totalVdx)
hold off


end

