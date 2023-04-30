function circplot(c, r, pars)
    t = linspace(0, 1, 100); 
    plot(c(1)+r*sin(2*pi*t), c(2)+r*cos(2*pi*t), pars)