clear all;
# x_1 - 3x_2 -> min
# -x_1 + 2x_2 <= 6
# x_1 + x_2 <= 5
# x_1, x_2 >= 0
c = [1; -3]; 
A = [-1, 2]; 
D = [1, 1]; 
b = [6];
d = [5];

lb = zeros(2, 1); ub = [];
ctype = repmat("U", 1, 2); 
vartype = repmat("C", 1, 2);

fprintf("####### Dantzig-Wolfle with known extreme points ######\n");
#{
all_extreme_points = [0, 0;
                      0, 5;
                      5, 0];
#}
start_of_poly_const = length(ctype) - size(D, 1)
ctype_sub = substr(ctype, start_of_poly_const + 1, size(D, 1))
all_extreme_points = extreme_points(D, d, lb, ub, ctype_sub, vartype);
fprintf("All extreme points:\n");
for i=1:size(all_extreme_points, 1)
  printf("(%d, %d) ", all_extreme_points(i, 1), all_extreme_points(i, 2));
endfor
printf("\n");
# odaberemo neke ekstremne točk1
extreme_points = [0,0; 0,5];
c_rest_master = [];
A_rest_master = [];
for i=1:size(extreme_points, 1)
  c_rest_master(i, :) = c'*extreme_points(i,:)';
  A_rest_master(:, i) = A*extreme_points(i,:)';
endfor  
# dodamo uvjet lamda_1 + lambda_2 = 1, lambda_1, lambda_2 >= 0,
# pa u matricu A dodamo redak 1, 1
A_rest_master(size(A_rest_master, 1)+1, :) = [1,1];
b_rest_master = b;
b_rest_master(size(b_rest_master, 1)+1, 1) = 1;
A_rest_master, c_rest_master, b_rest_master
# Restricted master
# 0*lambda_1 - 15*lambda_2 -> min
# 0*lambda_1 + 10*lambda_2 <= 6
# lambda_1 + lambda_2 = 1
# lambda_1, lambda_2 >= 0
# lower bound subproblem, upper bound subproblem
lb_rm = zeros(size(extreme_points, 1), 1); ub_rm = []; 
ctype_rest_master(1,1) = "U"; ctype_rest_master(1, 2) = "S";
vartype_rest_master = repmat("C", 1, 2);

printf("Restricted master solution:\n");
# find_ex_obj_val -> objective function value u LP-u za traženje ekstremne
# tocke koju cemo dodati (tj. stupca)
optimal_cost = intmin; tol = -0.001; iteration = 1;

ctype_extreme = ctype(2:size(ctype, 1), :);
vartype_extreme = vartype;
lbs = zeros(size(extreme_points, 1), 0); ubs = []; 

find_ex_obj_val = intmin;
total_num_of_extremes = size(all_extreme_points, 1);
num_of_variables = size(extreme_points, 1);
while (find_ex_obj_val < tol && num_of_variables < total_num_of_extremes)
  [lambda_opt, fval, ~, extra] = glpk(c_rest_master, A_rest_master, b_rest_master, 
                                  lb_rm, ub_rm, ctype_rest_master, vartype_rest_master, 1);
  printf("%d. iteration of DW:\n", iteration);
  printf('Optimal Objective Function Value = %f\n', fval);
  printf('Optimal Solution = '); disp(lambda_opt');
  printf('Dual:\n'); disp(extra.lambda);

  pi = extra.lambda(1,1);
  alpha = extra.lambda(2,1);
  pi, alpha
  
  printf("Finding next extreme point:\n");
  printf("Objective function coefficients of subproblem:\n");
  [extreme_point, find_ex_obj_val] = glpk(c - A'*pi, D, d, lbs, ubs, 
                                    ctype_extreme, vartype_extreme, 1);
                                      
  printf("Extreme point to be added: "); disp(extreme_point');
  printf("Object value of subproblem to find ext. point to be added: %f\n", find_ex_obj_val);
  optimal_cost = find_ex_obj_val - alpha;
  
  if(optimal_cost < tol)
    extreme_points = [extreme_points; extreme_point'];
    new_column = [A*extreme_point; 1];
    A_rest_master = [A_rest_master, new_column];
    c_rest_master = [c_rest_master; c'*extreme_point];
    num_of_variables = num_of_variables + 1;
    vartype_rest_master = repmat("C", 1, num_of_variables);
    lb_rm =  zeros(num_of_variables, 1);
  endif

  iteration = iteration + 1;
end

printf("\nAfter generating columns:\n");
[lambda_opt, fval, ~, extra] = glpk(c_rest_master, A_rest_master, b_rest_master, 
                                  lb_rm, ub_rm, ctype_rest_master, vartype_rest_master, 1);
printf("Optimal Objective Function Value = %f\n", fval);
printf("Optimal Solution - lambdas: "); disp(lambda_opt');
printf("\nOptimal solution of original master problem: ");
x_optimal = lambda_opt'*extreme_points;
disp(x_optimal);
f_min = x_optimal(1,1) - 3 * x_optimal(1, 2);
printf("\nObjective function value of original master problem: %f\n", f_min);