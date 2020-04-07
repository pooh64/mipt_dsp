param = 1200:1:1625;
[alpha, D, beta2, gamma] = arrayfun(@smf, param);
plot(param, D);

[alpha, D, beta2, gamma] = smf(1561.83);
printf("D(1561.83) = %g\n", D);
[alpha, D, beta2, gamma] = smf(1530.33);
printf("D(1530.33) = %g\n", D);