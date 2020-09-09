# rješavamo problem:
# c^{T} x -> min ... c^{T}x je funkcija cilja (objective function)
# Ax = b, x >= 0 ... uvjeti jednakosti

c = [2; -1; 1; 0; 0; 0]; 
A = [3,  1,  1,  1,  0,  0; 
     1, -1,  2,  0,  1,  0; 
     1,  1, -1,  0,  0,  1]; 
b = [6; 1; 2]; % desna strana uvjeta jednakosti

# lb = array of lower bounds (u ovom slučaju su to 0 jer 0 <= x) x >= lb
# ub = array of upper bounds (u ovom slučaju nema gornjih granica) x <= ub
lb = zeros(6, 1); ub = [];

# ctype = polje koje sadrži oznake za vrste uvjeta
# ctype: "F" - free, ignoriran uvjet; 
#        "U" - nejednakost s gornjom oznakom A_i*x <= b_i
#        "S" -  jednakost A_i*x = b_i
#        "L" - nejednakost s gornjom oznakom A_i * x >= b_i  
#        "D" - nejednakost s obje granice
ctype = repmat("S", 1, 3); 
# vartype = polje vrsta varijabli: "C" - kontinuirana varijabla;
#                                   "I" - cjelobrojno rješenje  


vartype = repmat("C", 1, 6);
# posljednji parametar je sense: 1 = min, -1 = max, 1 default
[x, fval, errnum] = glpk(c, A, b, lb, ub, ctype, vartype, -1);
fprintf('Solution Status = %d\n', errnum);
fprintf('Optimal Objective Function Value = %f\n', fval);
fprintf('Optimal Solution = '); disp(x')