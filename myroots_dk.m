function r = myroots_dk(a, tol)
    % Find the roots of a polynomial given its coefficients
    % a is a vector of coefficients in descending order of degree.
    % tol is a tolerance value example: tol = 1e-6;
    % Initialize the roots as equally spaced points on the unit circle

    % Initialize the roots as equally spaced points on the unit circle
    n = length(a) - 1;
    w = exp(2i*pi*(0:n-1)/n)'; % equally spaced points on the unit circle
    r = 0.5*sum(a(2:end))*w; % radius scaler
    
    % Iterate the root-finding algorithm until convergence
    while max(abs(mypolyval(a, r))) > tol
        for k = 1:n
            r(k) = r(k) - mypolyval(a, r(k))./prod(r(k) - r([1:k-1 k+1:n])); %  update number k value 
        end
    end
        % Iterate the root-finding algorithm until convergence
%     while max(abs(mypolyval(a, r))) > tol
%         for k = 1:n
%             p = poly(r([1:k-1 k+1:n])); % compute polynomial using all roots except k-th
%             r(k) = r(k) - mypolyval(a, r(k))/mypolyval(p, r(k)); %  update number k value 
%         end
%     end
end
