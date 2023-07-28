% Decode reaction information : reactions, complexes, species, parameters, signs

% Type of reactions:
% (1) forward 
%     A --> B or A --> B+C or A --> B+C+D ==> kf*A  
%     A+B --> C ==> kf*A*B 
% (2) bound A+B <--> C or C <--> A+B ==> kf*A*B-kr*C or kf*C-kr*A*B 
% (*) particular case: some species are named with the symbol '*' as [A*]
% NB: does not exist A+B <--> C+D
% NB: A SINGLE REACTION IS A SINGLE ARROW! 
%     Therefore the number of reactions equals the number of parameters.

% Degradation reaction:
% A --> null ==> the null species can be ignored, 
%                but the null complex need to be added
% NB: null --> A is not present; hence, there are no constant parts in ODEs

% INPUTs: 'Reaction_arrow' reactions with arrows (as graph)
%         'ReactionFlux_rate' algebraic expressions of the reactions
%         'species_names' name of the species
%         'c_names' name of the parameters

% OUTPUTs: 'ReactionFlux_decoded' ReactionFlux splitted in signs, coefficients and species  
%          'ReactionFlux_decoded_ind' ReactionFlux splitted in signs, coefficients and species as indexes
%          'species_products' list of all products of species (quadratic)
%          'species_products_ind' related indexes of species in products
%          'complexes' list of all complexes (nodes of the graph)
%          'complexes_ind' related indexes of species in complexes
%          'reactions' list of all reactions as single arrow, related flux, related complexes
%          'reactions_ind' related indexes of coefficients (=1), species (<=2), complexes (=2)
%          'nr_sigle_arrow' and 'nr_double_arrow' number of non-reversible/reversible reactions

function [ReactionFlux_decoded,ReactionFlux_decoded_ind,...
    species_products,species_products_ind,...
    complexes,complexes_ind,reactions,reactions_ind,...
    nr_single_arrow,nr_double_arrow, reactions2flux_rates] = ...
    decode_reactions(Reaction_arrow,ReactionFlux_rate,species_names,c_names)

% dimensions
n_fluxes = length(ReactionFlux_rate);
n_reactions = length(c_names);

ReactionFlux_decoded = cell(size(ReactionFlux_rate));
ReactionFlux_decoded_ind = cell(size(ReactionFlux_rate));
species_products = cell(n_fluxes,2); % #products <= #fluxes
species_products_ind = zeros(n_fluxes,2); % each product is composed by at most two species
complexes = cell(2*n_fluxes,1); % #complexes <= #fluxes 
complexes_ind = zeros(2*n_fluxes,3); % each complex is composed by at most three species
reactions = cell(n_reactions,6); % #reactions(->) [edges] = #parameters 
reactions_ind = zeros(n_reactions,5); % each reaction flux is composed by 1 cofficient, at most two species
                                      % each reaction is composed by at most two complexes
% index to store reactions as single arrows
ind_reaction = 0;
% count reversible and irreversible reactions
nr_single_arrow = 0;
nr_double_arrow = 0;

for i = 1:n_fluxes

    %% DECODE species, parameters, signs information from ReactionFlux_rate
    
    flux = ReactionFlux_rate{i}; % single ReactionFlux
    if i==516
        silvia=1;
    end
    minus = strfind(flux,'-'); % if minus isempty than (1), else (2)
    % NB: the 'minus' variable can be also ' - ' (check missing!)
    
    star = strfind(flux,'*'); % '*' delimits products
    
    % case (*) ------------------------------------------------------------
    % remove '*' inside '[ ]'   
    brk_lx = strfind(flux,'[');
    brk_rx = strfind(flux,']'); 
    for j=1:length(brk_lx)  
        star( (star>brk_lx(j) & star<brk_rx(j)) ) = [];
    end   
    %----------------------------------------------------------------------
       
    if isempty(minus) % case (1) ------------------------------------------
        
        nr_single_arrow = nr_single_arrow+1; % irreversible reaction
        
        coeff_lx = flux(1:star(1)-1); % first parameter       
        coeff_rx = []; % second parameter (none)

        spec1_rx = []; % first species to the rigth of '-' (none)
        spec2_rx = []; % second species to the rigth of '-' (none)
        
        % is there is no minus, star has length 1 if k*c, length 2 if k*c1*c2
        if length(star) == 1 % if there is a single species to the left of '-'
            
            spec1_lx = flux(star(1)+1:end); % first species left
            spec2_lx = []; % none second species left  
            
        elseif length(star) == 2 % if there are two species to the left of '-'

            spec1_lx = flux(star(1)+1:star(2)-1); % first species left
            spec2_lx = flux(star(2)+1:end); % second species left  

        end
        
        
        %------------------------------------------------------------------
        
    else % case (2) -------------------------------------------------------

        nr_double_arrow = nr_double_arrow+1; % reversible reaction
        
        star_lx = star(star<minus); % to the left of '-'
        star_rx = star(star>minus); % to the rigth of '-'
        
        coeff_lx = flux(1:star_lx(1)-1); % first parameter 
        coeff_rx = flux(minus+1:star_rx(1)-1); % second parameter 

        if length(star_lx) == 1 % if there is a single species to the left of '-' 
            
            spec1_lx = flux(star_lx(1)+1:minus-1); % first species left
            spec2_lx = []; % none second species left 
            
        elseif length(star_lx) == 2 % if there are two species to the left of '-'
            
            spec1_lx = flux(star_lx(1)+1:star_lx(2)-1); % first species left
            spec2_lx = flux(star_lx(2)+1:minus-1); % second species left

        end
        
        if length(star_rx) == 1 % if there is a single species to the rigth of '-'  
            
            spec1_rx = flux(star_rx(1)+1:end); % first species right
            spec2_rx = []; % none second species right
            
        elseif length(star_rx) == 2 % if there are two species to the right of '-'
            
            spec1_rx = flux(star_rx(1)+1:star_rx(2)-1); % first species right
            spec2_rx = flux(star_rx(2)+1:end); % second species right 
    
        end
                      
    end
    %----------------------------------------------------------------------    
    
    % case (*) ------------------------------------------------------------    
    spec1_lx(spec1_lx == '[' | spec1_lx == ']') = [];
    spec2_lx(spec2_lx == '[' | spec2_lx == ']') = [];
    spec1_rx(spec1_rx == '[' | spec1_rx == ']') = [];
    spec2_rx(spec2_rx == '[' | spec2_rx == ']') = []; 
    %----------------------------------------------------------------------
    
    % store the products of species
    if isempty(minus)
        if length(star) == 2
            species_products{i,1} = spec1_lx;
            species_products{i,2} = spec2_lx;          
        end
    else
        if length(star_lx) == 2
            species_products{i,1} = spec1_lx;
            species_products{i,2} = spec2_lx;
        end
        if length(star_rx) == 2  
            species_products{i,1} = spec1_rx;
            species_products{i,2} = spec2_rx;
        end
    end
    
    %----------------------------------------------------------------------
    % ReactionFluxes_decoded{i} = {'+',coeff_lx,spec1_lx,spec2_lx,'-',coeff_rx,spec1_rx,spec2_rx};
    ReactionFlux_decoded{i,1} = '+';
    ReactionFlux_decoded{i,2} = coeff_lx;
    ReactionFlux_decoded{i,3} = spec1_lx;
    ReactionFlux_decoded{i,4} = spec2_lx;   
    if isempty(minus) % case (1)
        ReactionFlux_decoded{i,5} = [];
    else % case (2)
        ReactionFlux_decoded{i,5} = '-';
    end
    ReactionFlux_decoded{i,6} = coeff_rx;
    ReactionFlux_decoded{i,7} = spec1_rx;
    ReactionFlux_decoded{i,8} = spec2_rx;
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------    
    % index of 'coeff1' in 'c_names' 
    coeff_lx_ind = find(strcmp(c_names,coeff_lx));
    % index of 'coeff2' in 'c_names' 
    coeff_rx_ind = find(strcmp(c_names,coeff_rx));
    
    % index of 'spec1_lx' in 'species_names'
    spec1_lx_ind = find(strcmp(species_names,spec1_lx));
    % index of 'spec2_lx' in 'species_names'
    spec2_lx_ind = find(strcmp(species_names,spec2_lx));
    
    % index of 'spec1_rx' in 'species_names'
    spec1_rx_ind = find(strcmp(species_names,spec1_rx));
    % index of 'spec2_rx' in 'species_names'
    spec2_rx_ind = find(strcmp(species_names,spec2_rx));
    
    % ReactionFluxes_ind{i} = {'+',coeff_lx_ind,spec1_lx_ind,spec2_lx_ind,'-',coeff_rx_ind,spec1_rx_ind,spec2_rx_ind};
    ReactionFlux_decoded_ind{i,1} = '+';
    ReactionFlux_decoded_ind{i,2} = coeff_lx_ind;
    ReactionFlux_decoded_ind{i,3} = spec1_lx_ind;
    ReactionFlux_decoded_ind{i,4} = spec2_lx_ind;
    if isempty(minus) % case (1)
        ReactionFlux_decoded_ind{i,5} = [];
    else % case (2)
        ReactionFlux_decoded_ind{i,5} = '-';
    end
    ReactionFlux_decoded_ind{i,6} = coeff_rx_ind;
    ReactionFlux_decoded_ind{i,7} = spec1_rx_ind;
    ReactionFlux_decoded_ind{i,8} = spec2_rx_ind;
    
    % index of 'species_products{i}' in 'species_names'
    if isempty(species_products{i,1}) || isempty(species_products{i,2})
        species_products_ind(i,1) = 0;
        species_products_ind(i,2) = 0;
    else
        species_products_ind(i,1) = find(strcmp(species_names,species_products{i,1}));
        species_products_ind(i,2) = find(strcmp(species_names,species_products{i,2}));
    end
    %----------------------------------------------------------------------
    
    %% DECODE reactions and complexes information from Reaction_arrow %%%%%
    
    reaction = Reaction_arrow{i};
  
    arrow_bound = strfind(reaction,' <-> ');
    arrow_forward = strfind(reaction,' -> '); % OBs: included degradation reactions
    
    % determine complexes in reaction, reactions as single arrows
    if ~isempty(arrow_forward) % case (1) ---------------------------------
        
        % COMPLEXES
        complex_lx = reaction(1:arrow_forward-1); % to the left of ' -> '
        complex_rx = reaction(arrow_forward+length(' -> '):end); % to the right of ' -> ' 
        
        % REACTIONS
        ind_reaction = ind_reaction+1;
        % store reaction arrow
        reactions{ind_reaction,1} = reaction;
        % store coefficients and species
        reactions{ind_reaction,2} = coeff_lx;
        reactions{ind_reaction,3} = spec1_lx;
        reactions{ind_reaction,4} = spec2_lx;
        % store indexes of coefficients and species
        reactions_ind(ind_reaction,1) = coeff_lx_ind;
        reactions_ind(ind_reaction,2) = spec1_lx_ind;
        if ~isempty(spec2_lx_ind)
            reactions_ind(ind_reaction,3) = spec2_lx_ind;
        end
        % store complexes
        reactions{ind_reaction,5} = complex_lx;
        reactions{ind_reaction,6} = complex_rx;
        
        % store association between forward and backward reactions
        reactions2flux_rates(ind_reaction, 1) = i;
        
    end
    if ~isempty(arrow_bound) % case (2) -----------------------------------   
                             % NB: need to check for constant species
        
        % COMPLEXES
        complex_lx = reaction(1:arrow_bound-1); % to the left of ' <-> '
        complex_rx = reaction(arrow_bound+length(' <-> '):end); % to the right of ' <-> '
        
        % REACTIONS
        ind_reaction = ind_reaction+2;
        % store reaction arrow
        reactions{ind_reaction-1,1} = strcat(reaction(1:arrow_bound-1),{' -> '},reaction(arrow_bound+length(' <-> '):end));
        reactions{ind_reaction,1} = strcat(reaction(arrow_bound+length(' <-> '):end),{' -> '},reaction(1:arrow_bound-1)); 
        % store coefficients and species
        reactions{ind_reaction-1,2} = coeff_lx;
        reactions{ind_reaction-1,3} = spec1_lx;
        reactions{ind_reaction-1,4} = spec2_lx;
        reactions{ind_reaction,2} = coeff_rx;
        reactions{ind_reaction,3} = spec1_rx;
        reactions{ind_reaction,4} = spec2_rx;
        
        % store indexes of coefficients and species
        
        %disp(ind_reaction);
        reactions_ind(ind_reaction-1,1) = coeff_lx_ind;
        reactions_ind(ind_reaction-1,2) = spec1_lx_ind;
        if ~isempty(spec2_lx_ind)
            reactions_ind(ind_reaction-1,3) = spec2_lx_ind;
        end
        reactions_ind(ind_reaction,1) = coeff_rx_ind;
        reactions_ind(ind_reaction,2) = spec1_rx_ind;
        if ~isempty(spec2_rx_ind)
            reactions_ind(ind_reaction,3) = spec2_rx_ind;
        end 
        % store complexes
        reactions{ind_reaction-1,5} = complex_lx;
        reactions{ind_reaction-1,6} = complex_rx;
        reactions{ind_reaction,5} = complex_rx;
        reactions{ind_reaction,6} = complex_lx;
        
        % store association between forward and backward reactions
        reactions2flux_rates(ind_reaction-1, 1) = i;
        reactions2flux_rates(ind_reaction, 1) = i;
        
    end
    
    % store complexes 
    complexes(2*(i-1)+1:2*i) = {complex_lx;complex_rx};
    % check if the two complexes in reaction i>1 already exist
    if i>1

        exist_complex_lx = find(strcmp(complexes(1:2*(i-1)),complex_lx),1);
        exist_complex_rx = find(strcmp(complexes(1:2*(i-1)),complex_rx),1);
        if ~isempty( exist_complex_lx )   
            complexes(2*(i-1)+1) = {0}; 
        end
        if ~isempty( exist_complex_rx )    
            complexes(2*i) = {0};       
        end
        
    end
    
    % determine indexes of species inside complexes         
    if ~isequal(complexes{2*(i-1)+1},0)
        
        % case (*) ------------------------------------------------------------
        star_complex = strfind(complex_lx,'*'); % '*' inside species names
        if ~isempty(star_complex)
            % remove '[ ]'
            brk_lx = strfind(complex_lx,'[');
            brk_rx = strfind(complex_lx,']');
            complex_lx([brk_lx,brk_rx]) = [];
        end
        %----------------------------------------------------------------------
        
        plus_lx = strfind(complex_lx,' + ');
        % OBs: the left complexes are composed by at most two species (one '+')
        
        if isempty(plus_lx)
            complexes_ind(2*(i-1)+1,1) = find(strcmp(species_names,complex_lx));
            if (i==511)
                silvia=0;
            end
        else
            complexes_ind(2*(i-1)+1,1) = find(strcmp(species_names,complex_lx(1:plus_lx-1)));
            complexes_ind(2*(i-1)+1,2) = find(strcmp(species_names,complex_lx(plus_lx+length(' + '):end)));
        end
        
        
    end
    if ~isequal(complexes{2*i},0)
        % OBs: degradation reaction 'A --> null' implies a null complex at
        % right-hand side of an irreversible reaction, composed of just null  
        
        % case (*) ------------------------------------------------------------
        star_complex = strfind(complex_rx,'*'); % '*' inside species names
        if ~isempty(star_complex)
            % remove '[ ]'
            brk_lx = strfind(complex_rx,'[');
            brk_rx = strfind(complex_rx,']');
            complex_rx([brk_lx,brk_rx]) = [];
        end
        %----------------------------------------------------------------------
       
        plus_rx = strfind(complex_rx,' + ');
        % OBs: the right complexes are composed by at most three species (two '+')
        
        if isempty(plus_rx)
            if ~isempty(find(strcmp(species_names,complex_rx),1))
                complexes_ind(2*i,1) = find(strcmp(species_names,complex_rx));
            else 
                complexes_ind(2*i,1) = -1; % the null complex is decoded as -1
            end                
        elseif length(plus_rx) == 1
            complexes_ind(2*i,1) = find(strcmp(species_names,complex_rx(1:plus_rx-1)));
            complexes_ind(2*i,2) = find(strcmp(species_names,complex_rx(plus_rx+length(' + '):end)));
        elseif length(plus_rx) == 2
            complexes_ind(2*i,1) = find(strcmp(species_names,complex_rx(1:plus_rx(1)-1)));
            complexes_ind(2*i,2) = find(strcmp(species_names,complex_rx(plus_rx(1)+length(' + '):plus_rx(2)-1)));
            complexes_ind(2*i,3) = find(strcmp(species_names,complex_rx(plus_rx(2)+length(' + '):end)));
        end
        
    end
      
end

% remove empty raw cells and zero raws
species_products( all(cellfun('isempty',species_products),2),: ) = [];
species_products_ind( all(~species_products_ind,2), : ) = [];
complexes(cellfun(@(x) isequal(x, 0), complexes)) = [];
complexes_ind( all(~complexes_ind,2), : ) = [];

%% Store indexes of complexes in reactions

for r = 1:n_reactions  
    r_complex = reactions(r,5:6);    
    reactions_ind(r,4) = find(strcmp(complexes,r_complex(1)));
    reactions_ind(r,5) = find(strcmp(complexes,r_complex(2)));
end

end