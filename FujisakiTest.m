%% Fujisaki model

p.t = -.5:.001:2.75;

% Baseline
p.F0b = 90;

% Phrase Command
p.Ap = [.5 .2 -.2]';
p.T0 = [-.1 1.1 2.4]';
p.alpha = 2.5;

% Accent Command
p.Aa = [.4 .25 .5 .1]';
p.beta = 20;
p.gamma = .9;
p.T1 = [.1 .5 1.25 2.2]';
p.T2 = [.25 1.1 2.1 2.3]';


close all
obj = fujisakiModel(p);
obj = obj.assembleContour();
f = obj.plotF0Contour();
set(gcf, 'Position', get(0, 'Screensize'))
exportgraphics(f,"FujisakiFigure.png")
close(f)