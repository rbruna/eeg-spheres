function gamma = mymcs_gamma_yao2003_err ( headmodel, order )

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

% Creates constants for the radii and their ratios.
r1  = headmodel.r (1);
r4  = headmodel.r (3);
r23 = headmodel.r (2) / headmodel.r (3);
r34 = headmodel.r (3) / headmodel.r (4);
r43 = headmodel.r (4) / headmodel.r (3);
r32 = headmodel.r (3) / headmodel.r (2);
r21 = headmodel.r (2) / headmodel.r (1);


% Calculates the auxiliary constants.
k     = orders ./ ( orders + 1 );
p     = r34 .^ ( 2 * orders + 1 ) .* ( 1 - c34 ) - ( 1 + k .* c34 );
q     = ( 1 - c34 ) - r34 .^ ( 2 * orders + 1 ) .* ( 1 + c34 ./ k );
s     = p ./ q;
alpha = k .* r23 .^ ( 2 * orders + 1 ) .* ( c23 - 1 ) + s .* ( k .* c23 + 1 );
beta  = r23 .^ ( 2 * orders + 1 ) .* ( k + c23 ) + s .* ( 1 + c23 );
gamma = alpha ./ beta;

% Calculates the constant G.
G1    = 1 ./ ( r1 .^ ( 2 * orders + 1 ) );
G2    = ( 1 + s ) ./ ( 1 + k .* r43 .^ ( 2 * orders + 1 ) );
G3    = ( 1 + gamma ) ./ ( 1 + s .* r32 .^ ( 2 * orders + 1 ) );
G4    = ( 2 * orders + 1 ) ./ ( orders .* ( 1 - c21 ) + gamma .* r21 .^ ( 2 * orders + 1 ) .* ( orders + c21 .* ( orders + 1 ) ) );
G     = G1 .* G2 .* G3 .* G4;

% Calculates the spherical harmonics filter W at the outer layer.
W4    = G .* r4 .^ ( 2 * orders + 1 ) .* ( 1 + k );

% Calculates gamma.
% W4 = ( 2 .* l + 1 ) .^ 4 ./ gamma ./ l;
gamma = 1 ./ ( orders .* W4 ./ ( 2 .* orders + 1 ) .^ 4 );
