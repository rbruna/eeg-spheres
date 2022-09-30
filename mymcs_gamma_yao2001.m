function gamma = mymcs_gamma_yao2001 ( headmodel, order )

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

% Makes sure that there is a fourth layer.
if numel ( headmodel.r ) < 4
    headmodel.r      (4) = headmodel.r    (3);
    headmodel.cond   (4) = headmodel.cond (3);
end

% Creates constants for the conductivity ratios.
c21 = headmodel.cond (2) / headmodel.cond (1);
c23 = headmodel.cond (2) / headmodel.cond (3);
c34 = headmodel.cond (3) / headmodel.cond (4);

% Creates constants for the radii.
r1  = headmodel.r (1);
r2  = headmodel.r (2);
r3  = headmodel.r (3);
r4  = headmodel.r (4);


% Calculates the auxiliary constants.
k     = orders ./ ( orders + 1 );
p     = r3 .^ ( 2 * orders + 1 ) .* r4 .^ -( 2 * orders + 1 ) .* ( 1 - c34 ) - ( 1 + k .* c34 );
q     = r3 .^ -( 2 * orders + 1 ) .* ( 1 - c34 ) - r4 .^ -( 2 * orders + 1 ) .* ( 1 + c34 ./ k );
s     = p ./ q;
alpha = r2 .^ ( 2 * orders + 1 ) .* ( k .* r2 .^ ( 2 * orders + 1 ) .* ( c23 - 1 ) + s .* ( k .* c23 + 1 ) );
% beta  = r2 .^ ( 2 * l + 1 ) .* ( k + c23 ) + s .* ( 1 + c23 );
% beta  = r2 .^ ( 2 * l + 1 ) .* ( k + c23 ) + s .* ( 1 - c23 ); % ( 1 + c23 ) -> ( 1 - c23 )
beta  = r2 .^ ( 2 * orders + 1 ) .* ( k + c23 ) - s .* ( 1 - c23 ); % ( 1 - c23 ) -> ( c23 - 1 )
gamma = alpha ./ beta;
xhi   = k .* r4 .^ ( 2 * orders + 1 );

% Calculates the constants A, G.
An    = ( 2 .* orders + 1 ) .* ( 1 + gamma .* r1 .^ -( 2 * orders + 1 ) );
Ad1   = orders .* r1 .^ ( 2 * orders + 1 ) .* ( 1 - c21 );
Ad2   = gamma .* ( orders + c21 .* ( orders + 1 ) );
A     = An ./ ( Ad1 + Ad2 ) - r1 .^ -( 2 * orders + 1 );

G1n   = 1 + s .* r3 .^ -( 2 * orders + 1 );
G1d   = 1 + k .* r4 .^ ( 2 * orders + 1 ) .* r3 .^ -( 2 * orders + 1 );
G1    = G1n ./ G1d;
G2n   = 1 + gamma .* r2 .^ -( 2 * orders + 1 );
G2d   = 1 + s .* r2 .^ -( 2 * orders + 1 );
G2    = G2n ./ G2d;
G3n   = 2 * orders + 1;
G3d   = orders .* r1 .^ ( 2 * orders + 1 ) .* ( 1 - c21 ) + gamma .* ( orders + c21 .* ( orders + 1 ) );
G3    = G3n ./ G3d;
G     = G1 .* G2 .* G3;

% Calculates the spherical harmonics filter W at the outer layer.
W1    = A .* r1 .^ ( 2 * orders + 1 ) + 1; %#ok<NASGU>
W4    = G .* ( r4 .^ ( 2 * orders + 1 ) + xhi );

% Calculates gamma.
% W4 = ( 2 .* l + 1 ) .^ 4 ./ gamma ./ l;
gamma = 1 ./ ( orders .* W4 ./ ( 2 .* orders + 1 ) .^ 4 );
