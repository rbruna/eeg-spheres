function gamma = mymcs_gamma_bru_old ( headmodel, order )

% Returns the Gamma series expansion for EEG leadfields.
%
% Use as:
%   Gamma = mymcs_gamma ( headmodel, order );
%
% where:
%   headmodel   FieldTrip concentric spheres definition:
%       headmodel.r     Radius of the spheres.
%       headmodel.cond  Conductivity of each sphere.
%   order       Number of terms for the series (default 60). 
%
% The series expansion is valid for a sphere centered at origin.

% This implementation is adapted from:
%   Naess et al. 2017 Front. Hum. Neurosci. 2017.490.

% Based on FieldTrip 20160222 functions:
% * eeg_leadfield4_prepare by Robert Oostenveld

% Initializes the empty inputs.
if nargin < 2 || isempty ( order )
    order       = 60;
end

% Creates the vector of orders.
orders      = 1: order;

% Sorts the spheres from the smallest to the largest
[ ~, idx ]  = sort ( headmodel.r );
headmodel.r    = headmodel.r    ( idx );
headmodel.cond = headmodel.cond ( idx );


% Reserves memory for the X and A variables.
X = nan ( numel ( headmodel.r ), order );
A_ = nan ( numel ( headmodel.r ), order );


% Gets the value of X for the outer sphere.
X ( end, : ) = orders ./ ( orders + 1 );

% Goes from the second-to-outer to the inner sphere.
for sindex = numel ( headmodel.r ) - 1: -1: 1
    
    % Gets the first and previous values.
    Xn = X ( end, : );
    Xp = X ( sindex + 1, : );
    
    % Gets the current and previous radii.
    Rp = headmodel.r ( sindex + 1 );
    Ri = headmodel.r ( sindex );
    
    % Gets the current and previous conductivities.
    Sp = headmodel.cond ( sindex + 1 );
    Si = headmodel.cond ( sindex );
    
    % Gets the repeated term.
    num = Xn .* ( ( Ri / Rp ) .^ orders ) - Xp .* ( ( Rp / Ri ) .^ ( orders + 1 ) );
    den = ( ( Ri / Rp ) .^ orders ) + Xp .* ( ( Rp / Ri ) .^ ( orders + 1 ) );
    term = num ./ den;
    
    
    % Gets the value of X for the current sphere.
    num = Xn * ( Si / Sp ) - term;
    den = ( Si / Sp ) + term;
    
    X ( sindex, : ) = num ./ den;
end


% Gets the value of A_ for the inner sphere.
A_ ( 1, : ) = 1 ./ X ( 1, : );

% Goes from the second to the outer sphere.
for sindex = 2: numel ( headmodel.r )
    
    % Gets the previous value of A_.
    A_p = A_ ( sindex - 1, : );
    
    % Gets the previous and current value of X.
    Xp = X ( sindex - 1, : );
    Xi = X ( sindex, : );
    
    % Gets the current and previous radii.
    Rp = headmodel.r ( sindex - 1 );
    Ri = headmodel.r ( sindex );
    
    % Gets the value of A_ for the current sphere.
    num = A_p .* ( 1 + Xp );
    den = ( ( Rp / Ri ) .^ orders ) + Xi .* ( ( Ri / Rp ) .^ ( orders + 1 ) );
    
    A_ ( sindex, : ) = num ./ den;
end


% X1  = A1 ./ r1 .^ ( 2 * orders + 1 ); %#ok<NASGU>
% X2  = A2 ./ r1 .^ ( orders + 1 ) ./ r2 .^ ( orders ); %#ok<NASGU>
% X3  = A3 ./ r1 .^ ( orders + 1 ) ./ r3 .^ ( orders ); %#ok<NASGU>
% X4  = A4 ./ r1 .^ ( orders + 1 ) ./ r4 .^ ( orders );

A4 = A_ ( end, : );
U  = X ( end, : );
r1 = headmodel.r ( 1 );
r4 = headmodel.r ( end );


X4  = A4 ./ r1 .^ ( orders + 1 ) ./ r4 .^ ( orders );

% Calculates the spherical harmonics filter X at the outer layer.
X4  = X4 .* ( r4 .^ ( 2 * orders + 1 ) + U .* r4 .^ ( 2 * orders + 1 ) );

% Calculates gamma.
% X4 = ( 2 .* l + 1 ) .^ 4 ./ gamma ./ l;
gamma = 1 ./ ( orders .* X4 ./ ( 2 .* orders + 1 ) .^ 4 );
