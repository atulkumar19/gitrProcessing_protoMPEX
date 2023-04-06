function [x_LB,x_UB,f_LB,f_MID,f_UB,minus_percent,plus_percent] =  worstcase_example
x = [4     5     6; ...
     9     1     0; ...
     23     2     1];
errors =   [1          0.5         0.25; ...
            2          0.2         1.00; ...
            1          1.2         1.00];

[x_LB,x_UB,f_LB,f_MID,f_UB,minus_percent,plus_percent]  = worstcase(@(x) f_example_worstcase(x),x,errors);
    function y = f_example_worstcase(x)
        y = (x(1).^2 - 3*x(2) + 5*x(3))/sqrt(x(1) + x(2) + x(3) + 1);
    end
end
