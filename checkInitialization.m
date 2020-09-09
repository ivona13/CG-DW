function [initialized, f] = checkInitialization(points, c, A, b, lb, ub, 
                          ctype, vartype, sense)
  
  number_of_variables = size(points, 2);
  number_of_constraints = size(A, 1);
  number_of_lambdas = size(points, 1);
  
  #[lambda1, labmda2]
  A_lambda = zeros(number_of_constraints+1, 2);
  
  for cons_ind = 1:number_of_constraints+1
    for lambda = 1:number_of_lambdas
      if cons_ind > number_of_constraints
        A_lambda(cons_ind, lambda) = 1;
      else
        A_lambda(cons_ind, lambda) = A(cons_ind, :)*points(lambda,:)';
      endif
    endfor
  endfor
  
  c_lambda = zeros(1, number_of_lambdas);
  #c, points(1, :)
  for lambda = 1:number_of_lambdas
    c_lambda(1, lambda) = c'*points(lambda, :)';
  endfor
 
  # dodaj uvjet lambda1 + lambda2 = 1
  b_lambda = [b; 1];
  lb_lambda = zeros(number_of_lambdas, 1);
  ub_lambda = [];
  ctype_lambda = [ctype(1, 1:number_of_constraints), "S"];
  vartype_lambda = repmat("C", 1, number_of_lambdas);
 
 [lambda_opt, fval] = glpk(c_lambda, A_lambda, b_lambda, lb_lambda, 
                                  ub_lambda, ctype_lambda, vartype_lambda, sense);
  
  if isnan(fval)
    initialized = false;
  else
    initialized = true;
  endif  
  f = fval;
endfunction 