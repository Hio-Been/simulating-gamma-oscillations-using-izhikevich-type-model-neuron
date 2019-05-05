function kernel_func = get_exp_decreasing_kernel( distance, expF )

dd = exp(linspace(expF, 0, max(distance))); dd=dd/max(dd);
[xData, yData] = prepareCurveData( distance, dd );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.StartPoint = [ 1, .08, 0 ];
[kernel_func, gof] = fit( xData, yData, ...
    fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' ), opts );
