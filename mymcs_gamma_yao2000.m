function gamma = mymcs_gamma_yao2000 ( headmodel, order )

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
%   Yao 2000 Clin. Neurophisiol. 111: 81-92.
%   Yao 2001 Phys. Med. Biol. 46: 3177-3189.
%   Yao 2003 Phys Med. Biol. 48: 1997-2011.

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

% Checks that there are only three layers.
if numel ( headmodel.r ) > 3
    error ( 'This method only works for three leyers.' )
end

% Checks that the first and third conductivities are equal.
if headmodel.cond (1) ~= headmodel.cond (3)
    error ( 'This method only works for simmetrical conductivities.' )
end

% Creates a constant for the conductivity ratio.
c23 = headmodel.cond (2) / headmodel.cond (3);

% Creates constants for the radii and their ratios.
r1  = headmodel.r (1);
r2  = headmodel.r (2);
r3  = headmodel.r (3);
r23 = headmodel.r (2) / headmodel.r (3);


% Calculates the auxiliary constants.
xhi   = orders ./ ( orders + 1 ) .* r3 .^ ( 2 * orders + 1 );
alpha = r23 .^ ( 2 * orders + 1 ) .* ( 1 - c23 ) - ( 1 + ( orders .* c23 ) ./ ( orders + 1 ) );
beta  = r23 .^ -( 2 * orders + 1 ) .* ( 1 - c23 ) - ( 1 + ( orders + 1 ) ./ orders .* c23 );
gamma = alpha ./ beta;

% Calculates the constants A, B, E.
An    = ( 2 * orders + 1 ) .* ( 1 + gamma .* r1 .^ -( 2 * orders + 1 ) );
Ad1   = orders .* ( 1 - c23 ) .* r1 .^ ( 2 * orders + 1 );
Ad2   = gamma .* ( orders + c23 .* ( orders + 1 ) );
A     = An ./ ( Ad1 + Ad2 ) - r1 .^ -( 2 * orders + 1 );

Bn    = 2 * orders + 1;
Bd1   = orders .* ( 1 - c23 ) .* r1 .^ ( 2 * orders + 1 );
Bd2   = gamma .* ( orders + c23 .* ( orders + 1 ) );
B     = Bn ./ ( Bd1 + Bd2 );

En    = orders + 1;
Ed1   = ( orders + 1 ) .* r2 .^ ( 2 * orders + 1 );
Ed2   = orders .* r3 .^ ( 2 * orders + 1 );
E     = B .* ( r2 .^ ( 2 * orders + 1 ) + gamma ) .* En ./ ( Ed1 + Ed2 );

% Calculates the spherical harmonics filter K at the outer layer.
K1    = A .* r1 .^ ( 2 * orders + 1 ) + 1; %#ok<NASGU>
K2    = B .* ( r2 .^ ( 2 * orders + 1 ) + gamma ); %#ok<NASGU>
K3    = E .* ( r3 .^ ( 2 * orders + 1 ) + xhi );

% Calculates gamma.
% K3 = ( 2 .* l + 1 ) .^ 4 ./ gamma ./ l;
gamma = 1 ./ ( orders .* K3 ./ ( 2 .* orders + 1 ) .^ 4 );
