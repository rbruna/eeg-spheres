function gamma = mymcs_gamma_lutkenhoner ( headmodel, order )

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
%   Lutkenhoner, Habilschrift 1992. Möglichkeiten und Grenzen der neuromagnetischen Quellenanalyse.
% which again is taken from:
%   Cuffin & Cohen 1979 Electroencephalogr. Clin. Neurophysiol. 47:131-146.

% Based on FieldTrip 20160222 functions:
% * eeg_leadfield4_prepare by Robert Oostenveld

% Copyright (C) 2016, Ricardo Bruna

% Initializes the empty inputs.
if nargin < 2 || isempty ( order )
    order       = 60;
end

% Creates the vector of orders.
orders      = 1: order;

% sort the spheres from the smallest to the largest
[ ~, idx ]  = sort ( headmodel.r );
headmodel.r    = headmodel.r    ( idx );
headmodel.cond = headmodel.cond ( idx );

% Makes sure that there is a fourth layer.
if numel ( headmodel.r ) < 4
    headmodel.r      (4) = headmodel.r    (3);
    headmodel.cond   (4) = headmodel.cond (3);
end

% Creates constants for the conductivity ratios.
% Cuffin & Cohen 1979. Eq A2.
k1 = headmodel.cond (1) / headmodel.cond (2);
k2 = headmodel.cond (2) / headmodel.cond (3);
k3 = headmodel.cond (3) / headmodel.cond (4);

% Creates extra constants for the radii ratios.
a = headmodel.r (1) / headmodel.r (2);
b = headmodel.r (1) / headmodel.r (3);
c = headmodel.r (2) / headmodel.r (3);
d = headmodel.r (3) / headmodel.r (4);

% The original defintion can be rewritten as (FieldTrip):
gamma1_1 = ( orders * k1 + orders + 1) .* ( orders * k2 + orders + 1 );
gamma1_2 = orders .* ( orders + 1 ) * ( k1 - 1 ) * ( k2 - 1 ) .* a .^ ( 2 * orders + 1 );
gamma1   = gamma1_1 + gamma1_2;
gamma2_1 = orders * k3 + orders + 1;
gamma2_2 = ( orders + 1 ) * ( k3 - 1 ) .* d .^ ( 2 * orders + 1 );
gamma2   = gamma2_1 + gamma2_2;
gamma3_1 = ( k1 - 1 ) * ( ( orders + 1 ) * k2 + orders ) .* b .^ ( 2 * orders + 1 );
gamma3_2 = ( orders * k1 + orders + 1 ) * ( k2 - 1 ) .* c .^ ( 2 * orders + 1 );
gamma3   = gamma3_1 + gamma3_2;
gamma4   = orders + 1;
gamma5_1 = orders * ( k3 - 1 );
gamma5_2 = ( ( orders + 1 ) * k3 + orders ) .* d .^ ( 2 * orders + 1 );
gamma5   = gamma5_1 + gamma5_2;

gamma    = gamma1 .* gamma2 + gamma3 .* gamma4 .* gamma5;
