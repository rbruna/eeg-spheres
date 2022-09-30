function gamma = mymcs_gamma_naess_old ( headmodel, order )

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

% Makes sure that there is a fourth layer.
if numel ( headmodel.r ) < 4
    headmodel.r      (4) = headmodel.r    (3);
    headmodel.cond   (4) = headmodel.cond (3);
end

% Calculates the conductivity ratios.
c12 = headmodel.cond (1) / headmodel.cond (2);
c23 = headmodel.cond (2) / headmodel.cond (3);
c34 = headmodel.cond (3) / headmodel.cond (4);

% Calculates the radii ratios.
r1  = headmodel.r (1);
r2  = headmodel.r (2);
r3  = headmodel.r (3);
r4  = headmodel.r (4);
r12 = headmodel.r (1) / headmodel.r (2);
r23 = headmodel.r (2) / headmodel.r (3);
r34 = headmodel.r (3) / headmodel.r (4);
r43 = headmodel.r (4) / headmodel.r (3);
r32 = headmodel.r (3) / headmodel.r (2);
r21 = headmodel.r (2) / headmodel.r (1);


% Calculates the (re-defined) auxiliary constants U, V, Y and Z.
U   = orders ./ ( orders + 1 );

Vn  = U .* c34 - ( U .* r34 .^ orders - U .* r43 .^ ( orders + 1 ) ) ./ ( r34 .^ orders + U .* r43 .^ ( orders + 1 ) );
Vd  = c34 + ( U .* r34 .^ orders - U .* r43 .^ ( orders + 1 ) ) ./ ( r34 .^ orders + U .* r43 .^ ( orders + 1 ) );
V   = Vn ./ Vd;

Yn  = U .* c23 - ( U .* r23 .^ orders - V .* r32 .^ ( orders + 1 ) ) ./ ( r23 .^ orders + V .* r32 .^ ( orders + 1 ) );
Yd  = c23 + ( U .* r23 .^ orders - V .* r32 .^ ( orders + 1 ) ) ./ ( r23 .^ orders + V .* r32 .^ ( orders + 1 ) );
Y   = Yn ./ Yd;

Zn  = U .* r12 .^ orders - Y .* r21 .^ ( orders + 1 );
Zd  = r12 .^ orders + Y .* r21 .^ ( orders + 1 );
Z   = 1 ./ U .* Zn ./ Zd;

% Calculates the (re-defined) coefficients As.
A1n = ( orders + 1 ) ./ orders * c12 + Z;
A1d = c12 - Z;
A1  = A1n ./ A1d;

A2n = A1 + 1;
A2d = r12 .^ orders + r21 .^ ( orders + 1 ) .* Y;
A2  = A2n ./ A2d;

A3n = A2 .* ( 1 + Y );
A3d = r23 .^ orders + r32 .^ ( orders + 1 ) .* V;
A3  = A3n ./ A3d;

A4n = ( orders + 1 ) ./ orders .* A3 .* ( 1 + V );
A4d = ( orders + 1 ) ./ orders .* r34 .^ orders + r43 .^ ( orders + 1 );
A4  = A4n ./ A4d;

% X1  = A1 ./ r1 .^ ( 2 * orders + 1 ); %#ok<NASGU>
% X2  = A2 ./ r1 .^ ( orders + 1 ) ./ r2 .^ ( orders ); %#ok<NASGU>
% X3  = A3 ./ r1 .^ ( orders + 1 ) ./ r3 .^ ( orders ); %#ok<NASGU>
X4  = A4 ./ r1 .^ ( orders + 1 ) ./ r4 .^ ( orders );

% Calculates the spherical harmonics filter X at the outer layer.
X4  = X4 .* ( r4 .^ ( 2 * orders + 1 ) + U .* r4 .^ ( 2 * orders + 1 ) );

% Calculates gamma.
% X4 = ( 2 .* l + 1 ) .^ 4 ./ gamma ./ l;
gamma = 1 ./ ( orders .* X4 ./ ( 2 .* orders + 1 ) .^ 4 );
