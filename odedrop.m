function [dydx] = odedrop(s,Y,curv,capil)
    capil = 13.45;
    curv = 10;
    dydx = [2*curv+capil*Y(3)-sin(Y(1))/Y(2); -cos(Y(1)); sin(Y(1))];

end