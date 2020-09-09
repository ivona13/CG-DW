# -4x_1 - x_2 -6x_3 -> min
# 3x_1 + 2x_2 + 4x_3 = 17
# x_1 <= 2
# x_2 <= 2
# x_3 <= 2
# x_1, x_2, x_3 >= 1
c = [-4; -1; -6]; 
A = [3, 2, 4]; 
D = [1 0 0;
    0, 1, 0;
    0, 0, 1]; 
b = [17];
d = [2; 2; 2];

lb = ones(3, 1); ub = [];
ctype = ["S", "U", "U", "U"]; 
vartype = repmat("C", 1, 3);

printf("####### Dantzig-Wolfle with known extreme points ######\n");
all_extreme_points = [1, 1, 1;
                      1, 1, 2;
                      1, 2, 1;
                      2, 1, 1;
                      1, 2, 2;
                      2, 1, 2;
                      2, 2, 1;
                      2, 2, 2];
printf("All extreme points:\n");
for i=1:size(all_extreme_points, 1)
  printf("(%d, %d, %d) ", all_extreme_points(i, 1), all_extreme_points(i, 2),
                          all_extreme_points(i, 3));
endfor
printf("\n");
# odaberemo neke ekstremne točke
extreme_points = [2, 1, 2; 1, 1, 2];
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

# lower bound subproblem, upper bound subproblem
lb_rm = zeros(size(extreme_points, 1), 1); ub_rm = []; 
ctype_rest_master(1,1) = "S"; ctype_rest_master(1, 2) = "S";
vartype_rest_master = repmat("C", 1, 2);

printf("Restricted master solution:\n");
# find_ex_obj_val -> objective function value u LP-u za traženje ekstremne
# tocke koju cemo dodati (tj. stupca)
optimal_cost = intmin; tol = -0.001; iteration = 1;

ctype_extreme = ["U", "U", "U"];
vartype_extreme = vartype;
lbs = lb; ubs = []; 

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
  printf("Objective function coefficients of subproblem:\n"); disp(c-A'*pi);
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
  else
    break;
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