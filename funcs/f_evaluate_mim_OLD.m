%% TODO per velocizzare.
% Io calcolo tutto il sistema S*v e poi sostituisco le leggi di cons.
% Potrei già calcolare il sistema ridotto S_2*v.

function f_x = f_evaluate_mim(rate_constants, x, idx_basic_species, ...
                                Nl, rho, MIM)
                            
f_x = f_odefun_MIM(0, x, rate_constants, MIM, 'Sv');
f_x(idx_basic_species) = Nl*x - rho;

end