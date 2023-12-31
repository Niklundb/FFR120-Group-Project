Mt_Bruno = csvread("mt_bruno_elevation.csv");
x = linspace(-2500, 2500, 25);
y = linspace(-2500, 2500, 26);

[X, Y] = meshgrid(x,y);

[fitresult, gof] = createFit(X, Y, Mt_Bruno);

x1 = linspace(-2500, 2500, 100);
y1 = linspace(-2500, 2500, 100);

[X1, Y1] = meshgrid(x1, y1);

Z = fitresult(X1, Y1);
surf(X1, Y1, Z)


%%

nTimeSteps = 20;
nPoints = 1000;
windDirection = 30*pi/180;
dt = 5;
smalldt = 0.1;
windSpeed = 10;
fireSpeed = 0.2;
Ks = 1;
slope = 15*pi/180;

thetaList = linspace(0, 2 * pi, nPoints);
initFireBoundaryx = 5*cos(thetaList);
initFireBoundaryy = 5*sin(thetaList);


FireBoundaryx = zeros(nPoints, nTimeSteps+1);
FireBoundaryy = zeros(nPoints, nTimeSteps+1);

FireBoundaryx(:, 1) = initFireBoundaryx;
FireBoundaryy(:, 1) = initFireBoundaryy;
plot(FireBoundaryx(:, 1), FireBoundaryy(:, 1), 'b')
hold on

slopes = zeros(nPoints, nTimeSteps);

xi = zeros(1, 1000);

dy = gradient(FireBoundaryy(:, 1));
dx = gradient(FireBoundaryx(:, 1));

for i = 1:nTimeSteps

    for j = 1:nPoints

        x = FireBoundaryx(j, i);
        y = FireBoundaryy(j, i);
        h = fitresult(x, y);
        xPert = x+smalldt*dy(j);
        yPert = y-smalldt*dx(j);
        hpert = fitresult(xPert, yPert);
        dDist = norm([x, y] - [xPert, yPert]);
        dH = hpert - h;
        slopes(j, i) = atan(dH/dDist);
        if dx(j) > 0 && dy(j) > 0
            xi(j) = atan(-dx(j)/dy(j)) + 2*pi;
        elseif dx(j) < 0 && dy(j) > 0
            xi(j) = atan(-dx(j)/dy(j));
        elseif dx(j) > 0 && dy(j) < 0
            xi(j) = atan(-dx(j)/dy(j)) + pi;
        elseif dx(j) < 0 && dy(j) < 0
            xi(j) = atan(-dx(j)/dy(j)) + pi;
        end
        nu = xi(j);
        if nu < 0
            nu = nu +2*pi;
        end
        if (0 < nu) && (nu < pi/2) || (3*pi/2 < nu) && (nu < 2*pi)
            fireSpread = fireSpeed*Ks*exp(3.533*sign(slopes(j,i))*abs(tan(slopes(j,i))*cos(nu))^1.2*cos(nu))*exp(0.1783*windSpeed*cos(nu-windDirection));
        else
            fireSpread = fireSpeed*Ks*exp(-3.533*sign(slopes(j,i))*abs(tan(slopes(j,i))*cos(nu))^1.2*cos(nu))*exp(0.1783*windSpeed*cos(nu-windDirection));
        end
        FireBoundaryx(j, i+1) = FireBoundaryx(j, i) + cos(xi(j))*fireSpread*dt;
        FireBoundaryy(j, i+1) = FireBoundaryy(j, i) + sin(xi(j))*fireSpread*dt;
    end
    FireBoundaryx(end, i+1) = FireBoundaryx(1, i+1);
    FireBoundaryy(end, i+1) = FireBoundaryy(1, i+1);
    plot(FireBoundaryx(: ,i+1), FireBoundaryy(:, i+1), 'b')
end

