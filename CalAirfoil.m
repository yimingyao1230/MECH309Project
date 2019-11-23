function [dydx] = CalAirfoil(x)
    toc = 0.08; % thickness ratio
    yx = toc * (-2 * x^2 + 82*x -840); % circular-arc defination for airfoil
    dydx = toc * (-4 * x + 82); % derivitive of circular-arc defination for airfoil
end