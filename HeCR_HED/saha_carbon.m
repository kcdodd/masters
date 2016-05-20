function [z,ratios, z_avg, dev] = saha_carbon(ni, Te, iterations)

e=1.602E-19;
me=9.11E-31;
h=6.626E-34;
L3 = (h^2/(2*pi*me*Te))^(3/2);

%b = 2*e^3*sqrt(pi/Te);
Ecl=1.1E-28*(ni^(1/3));


ne = 3*ni;
pop = zeros(7,1);
oldpop = zeros(7,1);
ratios = zeros(7,1);
diff = zeros(7,1);

pop(1) = ni/7;
pop(2) = ni/7;
pop(3) = ni/7;
pop(4) = ni/7;
pop(5) = ni/7;
pop(6) = ni/7;
pop(7) = ni/7;

Ip = [0, 11.26, 24.38, 47.88, 64.49, 392.1, 490.0]*e;
g = [9, 6, 1, 2, 1, 2, 1];

for i= 1:iterations
	oldpop = pop;
	u = g(2)*2*exp(-(Ip(2)-Ip(1)-Ecl)/Te)/(g(1)*L3*ne);
	pop(2) = max(1,(ni-(pop(3)+pop(4)+pop(5)+pop(6)+pop(7)))*u/(1+u));
	ne = pop(2)+2*pop(3)+3*pop(4)+4*pop(5)+5*pop(6)+6*pop(7);
	u=g(3)*2*exp(-(Ip(3)-Ip(2)-Ecl)/Te)/(g(2)*L3*ne);
	pop(3) = max(1,(ni-(pop(1)+pop(4)+pop(5)+pop(6)+pop(7)))*u/(1+u));
	ne = pop(2)+2*pop(3)+3*pop(4)+4*pop(5)+5*pop(6)+6*pop(7);
	u=g(4)*2*exp(-(Ip(4)-Ip(3)-Ecl)/Te)/(g(3)*L3*ne);
	pop(4) = max(1,(ni-(pop(1)+pop(2)+pop(5)+pop(6)+pop(7)))*u/(1+u));
	ne = pop(2)+2*pop(3)+3*pop(4)+4*pop(5)+5*pop(6)+6*pop(7);
	u=g(5)*2*exp(-(Ip(5)-Ip(4)-Ecl)/Te)/(g(4)*L3*ne);
	pop(5) = max(1,(ni-(pop(1)+pop(2)+pop(3)+pop(6)+pop(7)))*u/(1+u));
	ne = pop(2)+2*pop(3)+3*pop(4)+4*pop(5)+5*pop(6)+6*pop(7);
	u = g(6)*2*exp(-(Ip(6)-Ip(5)-Ecl)/Te)/(g(5)*L3*ne);
	pop(6) = max(1,(ni-(pop(1)+pop(2)+pop(3)+pop(4)+pop(7)))*u/(1+u));
	ne = pop(2)+2*pop(3)+3*pop(4)+4*pop(5)+5*pop(6)+6*pop(7);
	u = g(7)*2*exp(-(Ip(7)-Ip(6)-Ecl)/Te)/(g(6)*L3*ne);
	pop(7) = max(1,(ni-(pop(1)+pop(2)+pop(3)+pop(4)+pop(5)))*u/(1+u));
	ne = pop(2)+2*pop(3)+3*pop(4)+4*pop(5)+5*pop(6)+6*pop(7);
	pop = (pop/sum(pop))*ni;

	diff = oldpop-pop;
	dev = max(abs(diff)/ni);
end

ratios = [pop(1)/ni, pop(2)/ni, pop(3)/ni, pop(4)/ni, pop(5)/ni, pop(6)/ni, pop(7)/ni];
z = 0:6;

z_avg=z*(ratios');
