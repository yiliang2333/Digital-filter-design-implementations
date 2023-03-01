function [C, D] = FilterSplit(A, B)
%将滤波器分解为2阶和1阶的串联

% rA = roots(A); A: numerator
% rB = roots(B); B: denominator

% lenA = length(rA);
% lenB = length(rB);

% gain = abs(A(1) / B(1)) ^ (1 / lenB);

[C, rC] = Split(A);
[D, rD] = Split(B);

lC = size(C);
lC = lC(1);
lD = size(D);
lD = lD(1);

if lD > lC
    for i=lC+1:lD
        C(i, :) = [1, 0, 0];
        rC(i) = 0;
    end
elseif lC > lD
    for i=lD+1:lC
        D(i, :) = [1, 0, 0];
        rD(i) = 0;
    end
    lD = lC;
end


%sort D, 离单位圆越近的极点优先级越高
dD = abs(1 - abs(rD));
[~, oD] = sort(dD);
rD = rD(oD);
D = D(oD, :);

UseC = rC * 0 + 1;

for i=1:lD
    if max(UseC) == 0
        break;
    end
    [~, m2] = max(UseC ./ (abs(rC - rD(i)) + 0.00001));
    order(i) = m2;
    UseC(m2) = 0;
end

C = C(order, :);
co = abs(A(1)) ^ (1/lD);
C = C * co;
if A(1) < 0
    C(1, :) = -C(1, :);
end
hold off;
% for i=1:lD
%     filter1(C(i, :), D(i, :));
%     hold on;
% end

end





function [y, r2] = Split(x)
r = roots(x);
len = length(r);

i = 1;
num = 0;
n = 1;

while i <= len
   if i < len && imag(r(i)) ~= 0
       %复数根
       output = real(conv([1, -r(i)], [1, -r(i + 1)]));
       
       if i == 1 && x(1) < 0
           output = - output;
       end
       num = num + 1;
       r2(num) = r(i);
       i = i + 2;       
       y(num, :) = output';
   else
       num = num + 1;
       y(num, :) = [1, -r(i), 0];
       r2(num) = r(i);
       i = i + 1;       
   end
   
   n = n + 1;
end

end