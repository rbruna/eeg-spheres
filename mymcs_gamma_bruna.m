function gamma = mymcs_gamma_bruna ( headmodel, order )

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

% This implementation is based on:
%   Cuffin & Cohen 1979 Electroencephalogr. Clin. Neurophysiol. 47:131-146.
%   Yao 2000 Clin. Neurophisiol. 111: 81-92.
%   Yao 2001 Phys. Med. Biol. 46: 3177-3189.
%   Yao 2003 Phys Med. Biol. 48: 1997-2011.
%   Naess et al. 2017 Front. Hum. Neurosci. 2017:490.
%   Bruña et al. 2022 biorXiv.

% Based on FieldTrip 20160222 functions:
% * eeg_leadfield4_prepare by Robert Oostenveld

% Initializes the empty inputs.
if nargin < 2 || isempty ( order )
    order       = 60;
end

% Creates the vector of orders.
orders      = 1: order;

% Sorts the spheres from the smallest to the largest.
[ ~, idx ]  = sort ( headmodel.r );
headmodel.r    = headmodel.r    ( idx );
headmodel.cond = headmodel.cond ( idx );




% Reserves memory for the X and A variables.
X = nan ( numel ( headmodel.r ), order );
A_ = nan ( numel ( headmodel.r ), order );

% Calculates the relative radii and conductivities.
Rr = headmodel.r (:) ./ headmodel.r (:)';
Sr = headmodel.cond (:) ./ headmodel.cond (:)';


% Gets the value of X for the outer sphere.
X ( end, : ) = orders ./ ( orders + 1 );

% Goes from the second-to-outer to the inner sphere.
for sindex = numel ( headmodel.r ) - 1: -1: 1
    
    % Gets the first and previous values.
    Xn = X ( end, : );
    Xp = X ( sindex + 1, : );
    
    % Gets the repeated term.
    num = Xn .* ( Rr ( sindex, sindex + 1 ) .^ orders ) - Xp .* ( Rr ( sindex + 1, sindex ) .^ ( orders + 1 ) );
    den = ( Rr ( sindex, sindex + 1 ) .^ orders ) + Xp .* ( Rr ( sindex + 1, sindex ) .^ ( orders + 1 ) );
    term = num ./ den;
    
    
    % Gets the value of X for the current sphere.
    num = Xn * Sr ( sindex, sindex + 1 ) - term;
    den = Sr ( sindex, sindex + 1 ) + term;
    
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
    
    % Gets the value of A_ for the current sphere.
    num = A_p .* ( 1 + Xp );
    den = ( Rr ( sindex - 1, sindex ) .^ orders ) + Xi .* ( Rr ( sindex, sindex - 1 ) .^ ( orders + 1 ) );
    
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
