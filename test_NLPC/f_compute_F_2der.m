function der_F = f_compute_F_2der(rate_constants, vm, Sm, Nl, idx_basic_species)
    [n_species, n_reactions] = size(Sm);
    der_F = zeros(n_species, n_species, n_species);
    % der_F(i,j,k) = 0 quando k>n-p, quindi vado a modificare solo se k<=n-p
    der_v = zeros(n_reactions, n_species, n_species);
    n_CLs = size(Nl,1);
    Sm(idx_basic_species, :) = [];
    for k=1:n_reactions
        for i=1:n_species
            for j=1:n_species
                if (vm(k,2) == vm(k,3)) && (i == j) && (i == vm(k,2))
                    der_v(k,i,j) = 2*rate_constants(vm(k,1));
                elseif (vm(k,2) ~= vm(k,3)) && (((i == vm(k,2)) && ...
                        (j == vm(k,3))) || ((j == vm(k,2)) && (i == vm(k,3))))
                    der_v(k,i,j) = 1*rate_constants(vm(k,1));
                end
            end
        end
    end
    for k=1:n_species-n_CLs
        for i=1:n_species
            for j=i:n_species
                der_F(k,i,j) = Sm(k,:) * der_v(:,i,j);
            end
        end
    end

end

