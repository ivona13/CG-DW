function [extreme_points] = extreme_points(A, b, lb, ub, ctype, vartype)
  # m uvjeta, n varijabli
  [m, n] = size(A);
  A_extended = A;
  b_extended = b;
  ctype_extended = ctype;
  
  #dodaj uvjete za x_i >= lb_i, x_i <= ub
  for var=1:n
    vector = zeros(1, n);
    vector(:, var) = 1;
    if size(lb, 1) >= var
      A_extended = [A_extended; vector];
      b_extended = [b_extended; lb(var,1)];
      ctype_extended = [ctype_extended, "L"];
    endif
    if size(ub, 1) >= var
      A_extended = [A_extended; vector];
      b_extended = [b_extended; ub(var, 1)];
      ctype_extended = [ctype_extended, "U"];
    endif
  endfor
  
  #A_extended, b_extended, ctype_extended
  
  potential_extreme_points = [];
  # biramo n uvjeta (s n varijabli) (A') i rješavamo sustav A'x=b'
  [m_, n_] = size(A_extended);
  permutations = nchoosek(1:m_, n_);
  number_of_combinations = size(permutations, 1);
  
  for i=1:number_of_combinations
    #generiraj A'
    A_i = zeros(n, n);
    b_i = zeros(n, 1);
    ctype_i = repmat("", 1, n);
    for j=1:n
      A_i(j, :) = A_extended(permutations(i, j), :);
      b_i(j, :) = b_extended(permutations(i, j), :);
      ctype_i(:, j) = ctype_extended(:, permutations(i, j));
    endfor
    if det(A_i > 0)
      x = A_i \ b_i;
      if ismember(x', potential_extreme_points, 'rows') < 1
        potential_extreme_points = [potential_extreme_points; x'];
      endif
    endif
  endfor
      
  extreme_points = [];    
  #provjeri jesu li tocke zaista ekstremne, tj. jesu li ujedno i dopustive 
  for ind=1:size(potential_extreme_points, 1)
    pot_extreme_point = potential_extreme_points(ind,:);
    if checkIfFeasible(pot_extreme_point, A, b, lb, ub, ctype, vartype) 
      extreme_points = [extreme_points; pot_extreme_point];
    endif
  endfor
  
endfunction