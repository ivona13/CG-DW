c = [2; -1; 1; 0; 0; 0]; 
A = [3,  1,  1,  1,  0,  0; 
     1, -1,  2,  0,  1,  0; 
     1,  1, -1,  0,  0,  1]; 
b = [6; 1; 2];

lb = zeros(6, 1); ub = [];
ctype = repmat("S", 1, 3); 
vartype = repmat("C", 1, 6);

[x, fval, errnum] = glpk(c, A, b, lb, ub, ctype, vartype, -1);
fprintf('Optimal Objective Function Value = %f\n', fval);
fprintf('Optimal Solution = '); disp(x')

# primarni problem: 2x_1 - x_2 + x_3 -> max
#                   3x_1 + x_2 + x_3 + x_4 = 6
#                   x_1 - x_2 + 2x_3 + x_5 = 1
#                   x_1 + x_2 - x_3 + x_6 = 2
#                   x >= 0

# dualni problem: 6y_1 + 1 y_2 + 2y_3 -> min
#                 3y_1 + y_2 + y_3 >= 2
#                  y_1 - y_2 + y_3 >= -1
#                  y_1 + 2y_3 - y_3 >= 1
#                  y_1 >= 0
#                  y_2 >= 0
#                  y_3 >= 0
#                  y_1, y_2, y_3 slobodni

# Dualni
lb_dual = []; ub_dual = [];

ctype_dual = repmat("L", 1, 6);
vartype_dual = repmat("C", 1, 3);

[y, fval] = glpk(b, A', c, lb_dual, ub_dual, ctype_dual, vartype_dual, 1);
fprintf('Optimal Objective Function Value for Dual Problem = %f\n', fval);
fprintf('Optimal Solution for Dual Problem = '); disp(y')