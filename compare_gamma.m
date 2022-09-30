clc
clear
close all

% Defines the paramaters for the simulation.
order = 100;
nrep  = 1e6;

% Generates the standard head model.
headmodel        = [];
headmodel.type   = 'concentricspheres';
headmodel.o      = [ 0.0000 0.0000 0.0000 ];
headmodel.r      = [ 0.8700 0.9200 1.0000 ];
headmodel.cond   = [ 1.0000 0.0125 1.0000 ];
headmodel.tissue = { 'brain' 'skull' 'scalp' };

% % Generates a random (symmetrical) head model.
% headmodel        = [];
% headmodel.type   = 'concentricspheres';
% headmodel.o      = [ 0.0000 0.0000 0.0000 ];
% headmodel.r      = [ 0.8700 0.9200 1.0000 ];
% headmodel.cond = rand ( 2, 1 );
% headmodel.cond = -sort ( -headmodel.cond );
% headmodel.cond (3) = headmodel.cond (1);
% headmodel.tissue = { 'brain' 'skull' 'scalp' };


% Calculates the gamma constants using the different methods.
tic
for c = 1: nrep
gamma_lutkenhoner = mymcs_gamma_lutkenhoner ( headmodel, order );
end
t = toc;
fprintf ( 1, 'Elapsed time for Lutkenhoner method is %.4f ms.\n', 1000 * t / nrep );

tic
for c = 1: nrep
gamma_yao2000     = mymcs_gamma_yao2000 ( headmodel, order );
end
t = toc;
fprintf ( 1, 'Elapsed time for Yao (2000) method is %.4f ms.\n', 1000 * t / nrep );

tic
for c = 1: nrep
gamma_yao2001     = mymcs_gamma_yao2001 ( headmodel, order );
end
t = toc;
fprintf ( 1, 'Elapsed time for Yao (2001) method is %.4f ms.\n', 1000 * t / nrep );

tic
for c = 1: nrep
gamma_yao2003     = mymcs_gamma_yao2003 ( headmodel, order );
end
t = toc;
fprintf ( 1, 'Elapsed time for Yao (2003) method is %.4f ms.\n', 1000 * t / nrep );

tic
for c = 1: nrep
gamma_naess       = mymcs_gamma_naess ( headmodel, order );
end
t = toc;
fprintf ( 1, 'Elapsed time for Naess (2017) method is %.4f ms.\n', 1000 * t / nrep );

tic
for c = 1: nrep
gamma_bruna       = mymcs_gamma_bruna ( headmodel, order );
end
t = toc;
fprintf ( 1, 'Elapsed time for Bruna (2022) method is %.4f ms.\n', 1000 * t / nrep );


% Creates the output figure.
figure ( 'Units', 'centimeters', 'Position', [  1.0  1.0 12.0  7.0 ] )
axes ( 'Units', 'centimeters', 'Position', [  0.8  0.5 10.8  6.0 ], 'NextPlot', 'add' )

% Plots the relative error.
plot ( ( gamma_lutkenhoner - gamma_lutkenhoner ) ./ ( gamma_lutkenhoner + gamma_lutkenhoner ), '.' )
plot ( ( gamma_yao2000 - gamma_lutkenhoner ) ./ ( gamma_yao2000 + gamma_lutkenhoner ), '.' )
plot ( ( gamma_yao2001 - gamma_lutkenhoner ) ./ ( gamma_yao2001 + gamma_lutkenhoner ), '+' )
plot ( ( gamma_yao2003 - gamma_lutkenhoner ) ./ ( gamma_yao2003 + gamma_lutkenhoner ), '*' )
plot ( ( gamma_naess - gamma_lutkenhoner ) ./ ( gamma_naess + gamma_lutkenhoner ), 'o' )
plot ( ( gamma_bruna - gamma_lutkenhoner ) ./ ( gamma_bruna + gamma_lutkenhoner ), 'v' )

% Plots the legend.
legend ( { 'Zero line' 'Yao 2000 vs. Lutkenhoner' 'Yao 2001 vs. Lutkenhoner' 'Yao 2003 vs. Lutkenhoner' 'Naess 2017 vs. Lutkenhoner' 'Bruna 2022 vs. Lutkenhoner' }, 'Box', 'off' )
ax = findall ( gcf, 'Type', 'legend' );
pos = get ( ax, 'Position' );
pos (2) = 1 - pos (4);
set ( ax, 'Position', pos )

% Sets the view span.
range = max ( abs ( ylim ) );
ylim ( range * [ -1.05 +1.60 ] )

% Saves the figure.
set ( gcf, 'Renderer', 'painters' )
print ( '-dpng', '-r300', 'Relative error.png' )
