function [isFeasible] = checkIfFeasible(point, A, b, lb, ub, ctype, vartype)
  number_of_variables = size(point, 2);
  number_of_constraints = size(A, 1);
  number_of_lb = size(lb,1);
  number_of_ub = size(ub, 1);
  
  #provjeri varijable
  for variable = 1:number_of_variables
      #provjeri lb;
      if variable <= number_of_lb && point(1, variable) < lb(variable, 1)
        isFeasible = false;
        return
      endif
      # provjeri ub;
      if variable <= number_of_ub && point(1, variable) > ub(variable, 1)
        isFeasible = false;
        return
      endif
    
      # provjeri je li tip varijable odgovarajuci 
      if vartype(1, variable) == 'I' && point(1, variable) ~= int64(point(1, variable))
        isFeasible = false;
        return
      endif
  endfor
    
  #provjeri uvjete
  for cons_ind = 1:number_of_constraints
    if(strcmp(ctype(1, cons_ind), "U") && A(cons_ind, :)*point' > b(cons_ind, 1))
      isFeasible = false;
      return
    endif
    
    if(strcmp(ctype(1, cons_ind), "S") && A(cons_ind, :)*point' ~= b(cons_ind, 1))
      isFeasible = false;
      return
    endif
    
    if(strcmp(ctype(1, cons_ind), "L") && A(cons_ind, :)*point' < b(cons_ind, 1))
      isFeasible = false;
      return
    endif
  endfor
  
  isFeasible = true;

endfunction