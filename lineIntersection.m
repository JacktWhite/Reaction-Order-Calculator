% Function to calculate the tangent line intersection
function [intercept] = lineIntersection(leftTangentLineEqn, rightTangentLineEqn)
    equation1_sym = str2sym(leftTangentLineEqn);
    equation2_sym = str2sym(rightTangentLineEqn);
    
    [sol_x, sol_y] = solve([equation1_sym, equation2_sym]);
    intercept = [sol_x, sol_y];
end