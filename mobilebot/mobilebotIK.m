function y = mobilebotIK(u)
%mobilebotIK Accept xyz vector u and return joint angle vector y
%Remember to convert from global frame to cart frame (subtract xyz of cart
%pos)

x = u(1);
y = u(2);
z = u(3) - 0.3; %remove L0 from calc

th1 = atan2(x,y);
R = sqrt(x^2+y^2);
r = sqrt(z^2+R^2);
gamma = acos(r);
alpha = acos(1-2*r^2);
beta = atan2(z,R);
th2 = -gamma - beta;
th3 = pi - alpha;

y = [ 0 0 0 th1 th2 th3]';
end

