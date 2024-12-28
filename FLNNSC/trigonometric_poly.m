function phi_x = trigonometric_poly(x, d)
%% construct Trigonometric polunomials mapping(second-order)
% T0(x)=1;T1(x)=x;T2(x)=2x^2-1

phi_x = [];
for jj = 1:d
    phi_x = [phi_x, [x(jj), sin(pi*x(jj)), cos(pi*x(jj)), sin(2*pi*x(jj)), cos(2*pi*x(jj))]];
end
end