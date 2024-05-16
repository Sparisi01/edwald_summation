data = readmatrix("../output/cpu_time_comparison.csv");

N = data(1,:)';
ed = data(2,:)';
cou = data(3,:);

fit = functionFit();
fit.datax = N;
fit.datay = ed;

fitter.par = 0;
fit.model = @(par, x) x^par(1);

fit.plotModelFit();
