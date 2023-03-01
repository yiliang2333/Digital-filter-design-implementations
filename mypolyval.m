function y = mypolyval(p, x)
    % Evaluates a polynomial at one or more input values.
    % p is a vector of coefficients in descending order of degree.
    % x is a scalar or vector of input values.
    % y is the value of the polynomial at the input value(s).
    
    % Determine the degree of the polynomial
    n = length(p) - 1;
    
    % Initialize the output vector
    y = zeros(size(x));
    
    % Compute the polynomial for each input value
    for i = 1:numel(x)
        % Initialize the polynomial value for this input value
        yval = p(1); % from the coefficient of  x^n to x^0
        
        % Compute the polynomial value using 秦九昭's method
        for j = 2:n+1
            yval = yval*x(i) + p(j);
        end
        
        % Store the polynomial value in the output vector
        y(i) = yval;
    end
    
end
