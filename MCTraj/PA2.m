function [] = PA2()
% Baldeep Kooner 101004107
% ELEC 4700 MonteCarlo PA (2) Assignment
% nTime = number of time steps
% Ex = Electric Field in the x-direction
Ex = 5*10^-20;
dt = 1; 
nTime = 100; 
Ps = 0.05; %probability of scattering
%P = @(dx) 1 - exp(-dx/lambda); 
%P = 1 - 1.2*rand();
%x0 = 0; 
%v0 = 0;
t = 1:nTime;%time
dt = 1; 
x = zeros(1, nTime); %distance
vx = zeros(1, nTime); %velocity
vdx = zeros(1, nTime);
vxsum = zeros(1, nTime);
for n = 2:nTime
    P = rand();
    if P <= 0.05
        vx(n) = 0;
    else
        vx(n) = vx(n-1) + (1.60217653e-19)*Ex*dt / (9.10938215e-31);
    end
     %x(n) = vx(n-1)*dt + ((1.60217653e-19)*Ex*dt^2) / (2*(9.10938215e-31));
     x(n) = x(n-1) + vx(n-1)*dt; 
     vxsum(n) = vxsum(n-1) + vx(n);
     vdx(n) = vxsum(n)/n;
end


    
%plot(t, vx)
%xlabel('time')
%ylabel('velocity')
%hold on
%plot(t, x)
    subplot(1, 2, 1)
    comet(t, x)
    xlabel('time')
    ylabel('distance')
    subplot(1, 2, 2)
    comet(t, vx)
    hold on
    comet(t, vdx)

    xlabel('time')
    ylabel('velocity')
    %title(vdx)
    hold off





end

