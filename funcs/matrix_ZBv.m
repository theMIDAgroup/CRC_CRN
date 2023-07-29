% Write ODE system as: d(X)/dt = S v = Z B v 
% n = #species  p = #parameters
% c = #complexes [nodes]  r = #reactions(->) [edges] = #parameters 
% where: X vector of all species [dim(X) = n]
%        S stoichiometric matrix [dim(S) = n x r]
%        v vector of fluxes [dim(v) = r x 1]
%        Z complex stoichiometric matrix [dim(Z) = n x c]
%        B incidence matrix [dim(B) = c x r]

% It holds:
% - Z B = S stoichiometric matrix

% Type of reactions:
% (1) forward 
%     A --> B or A --> B+C or A --> B+C+D ==> kf*A  
%     A+B --> C ==> kf*A*B 
% (2) bound A+B <--> C or C <--> A+B ==> kf*A*B-kr*C or kf*C-kr*A*B 
% NB: do not exist A+B <--> C+D
% NB: A SINGLE REACTION IS A SINGLE ARROW! 
%     Therefore the number of reactions equals the number of parameters.

% INPUTs: 'species' list of species
%         'complexes_ind' indexes of species composing the complexes (nodes of the graph)
%         'reactions_ind' indexes of coefficients and species composing each reaction (single arrow) 

% OUTPUTs: 'S' stoichiometric matrix
%          'v' matrix r x 3 of indexes of coefficients (one) and species (at most two) for each flux
%          'Z' matrix n x c in {0,1,2} s.t. Z(i,j)=z indicates how (z-times) species i appear in complex j
%          'B' matrix c x r in {0,1,-1} s.t. B(i,j)=b indicates if complex i is a reactant (-1) or product (+1) in reaction j 
%          'ind_one' new index to decode the linear products composed by just one species

function [S,v,Z,B,ind_one] = ...
    matrix_ZBv(species,complexes_ind,reactions_ind)

% system dimensions
n_reactions = size(reactions_ind,1);
n_complexes = size(complexes_ind,1);
n_ode_species = length(species);

% new index to decode the linear products composed by just one species
ind_one = length(species)+1;

% build v: it is just the first three columns of reactions_ind 
v = reactions_ind(:,1:3);
v(v==0) = ind_one;

% initializations
Z = zeros(n_ode_species,n_complexes);
B = zeros(n_complexes,n_reactions);

% build Z 
for c = 1:n_complexes    
    complex = complexes_ind(c,:);

    if complex(1)>0 % if the complex is not the null complex (decoded as -1)       
        % Z(i,:) indicates the complexes containing species of index i
        if complex(1)==complex(2)
            Z(complex(1),c) = 2; % if the complex is composed by the sum of the same two species
            % NB: the complex made of three species may
            % have two equal species in other positions (check missing!)
        else
            Z(complex(1),c) = 1;
            if complex(2)~=0 % if the complex is composed by two species
                Z(complex(2),c) = 1;
                if complex(3)~=0 % if the complex is composed by three species
                    Z(complex(3),c) = 1;
                end
            end
        end
    end
end

% build B
for r = 1:n_reactions
    % indexes of complexes in reaction r
    r_complex = reactions_ind(r,4:5);    
    B(r_complex(1),r) = -1; % reactant
    B(r_complex(2),r) = +1; % product   
end

% build S
S = Z*B;

end