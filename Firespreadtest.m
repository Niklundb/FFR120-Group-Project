clc
nTimeSteps = 5;
nPoints = 1000;
windDirection = 45*pi/180;
dt = 1;
windSpeed = 6;
fireSpeed = 0.2;
Ks = 1;
slope = 0*pi/180;
slopeDirection = 0;
initFireBoundaryx = 5*cos(thetaList);
initFireBoundaryy = 5*sin(thetaList);

FireBoundaryx = zeros(1000, nTimeSteps+1);
FireBoundaryy = zeros(1000, nTimeSteps+1);

FireBoundaryx(:, 1) = initFireBoundaryx;
FireBoundaryy(:, 1) = initFireBoundaryy;
plot(FireBoundaryx(:, 1), FireBoundaryy(:, 1))
hold on

xi = zeros(1, 1000);

for i = 1:nTimeSteps
    dy = gradient(FireBoundaryy(:, i));
    dx = gradient(FireBoundaryx(:, i));
    for j = 1:nPoints
        if dx(j) > 0 && dy(j) > 0
            xi(j) = atan(-dx(j)/dy(j)) + 2*pi;
        elseif dx(j) < 0 && dy(j) > 0
            xi(j) = atan(-dx(j)/dy(j));
        elseif dx(j) > 0 && dy(j) < 0
            xi(j) = atan(-dx(j)/dy(j)) + pi;
        elseif dx(j) < 0 && dy(j) < 0
            xi(j) = atan(-dx(j)/dy(j)) + pi;
        end
        nu = xi(j) - slopeDirection;
        if nu < 0
            nu = nu +2*pi;
        end
        if (0 < nu) && (nu < pi/2) || (3*pi/2 < nu) && (nu < 2*pi)
            fireSpread = fireSpeed*Ks*exp(3.533*(tan(slope)*cos(nu))^1.2*cos(nu))*exp(0.1783*windSpeed*cos(nu-windDirection));
        else
            fireSpread = fireSpeed*Ks*exp(-3.533*(tan(slope)*cos(nu))^1.2*cos(nu))*exp(0.1783*windSpeed*cos(nu-windDirection));
        end
        FireBoundaryx(j, i+1) = FireBoundaryx(j, i) + cos(xi(j))*fireSpread*dt;
        FireBoundaryy(j, i+1) = FireBoundaryy(j, i) + sin(xi(j))*fireSpread*dt;
    end
    plot(FireBoundaryx(:, i+1), FireBoundaryy(:, i+1))
end


