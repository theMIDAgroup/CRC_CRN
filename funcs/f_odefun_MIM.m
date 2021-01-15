% Write ODE system as: d(X)/dt = S v = Z B v 
% n = #species  p = #parameters
% c = #complexes [nodes]  r = #reactions(->) [edges] = #parameters 
% where: X vector of all species [dim(X) = n]
%        S stoichiometric matrix [dim(S) = n x r]
%        v vector of fluxes [dim(v) = r x 1]
%        Z complex stoichiometric matrix [dim(Z) = n x c]
%        B incidence matrix [dim(B) = c x r]

% 'odefun' is needed for Matlab ODE solver 'ode15s':
% dXdt = odefun(t,X), for a scalar t and a column vector X, 
% must return a column vector dXdt of data type single or double that 
% corresponds to dXdt = f(t,X);
% odefun must accept both input arguments, t and X, even if one of the 
% arguments is not used in the function.

% INPUTs: 't' time variable
%         'X' vector of ODE species
%         'coeff' vector of ODE system coefficients
%         'MIM' MIM struct containing all the needed MIM elements
%         'choise' string vector for the selection of the system type to solve ('Sv','ZBv')

% OUTPUTs: 'dXdt' vector of derivatives of species

function dXdt = f_odefun_MIM(t, X, coeff, MIM, choise)

% select the needed MIM elements
S = MIM.matrix.S;
% Z = MIM.matrix.Z;
% B = MIM.matrix.B;
v = MIM.matrix.v;
ind_one = MIM.matrix.ind_one;
ODEs_species_is_constant = MIM.species.is_constant;
% ODEs_species_is_constant =  MIM.species.is_constant(MIM.species.is_in_ODE);

% initialization
dXdt = zeros(size(X));

% add the 'one' element to X
X(ind_one) = 1;

% evaluate v in coeff and species values
coeff_eval = coeff(v(:,1));
species1_eval = X(v(:,2));
species2_eval = X(v(:,3));

% ODE system
switch choise
    case 'Sv'

        dXdt(~ODEs_species_is_constant) = S(~ODEs_species_is_constant,:) * (coeff_eval.*species1_eval.*species2_eval);

        if any(ODEs_species_is_constant)
            dXdt(ODEs_species_is_constant) = 0;
        end
        
    case strcmp(choise,'ZBv')

        dXdt(~ODEs_species_is_constant) = Z(~ODEs_species_is_constant,:) * B * (coeff_eval.*species1_eval.*species2_eval);
        if any(ODEs_species_is_constant)
            dXdt(ODEs_species_is_constant) = 0;
        end
        
    otherwise
        error('The choise for evluation of MIM ODEs should be ''Sv'' or ''ZBv''! ');

end

end