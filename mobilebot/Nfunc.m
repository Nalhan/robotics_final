function N = Nfunc(in1)
%NFUNC
%    N = NFUNC(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    10-May-2018 12:26:43

th5 = in1(5,:);
th6 = in1(6,:);
t2 = th5+th6;
t3 = cos(t2);
N = [0.0;0.0;1.566e3;0.0;t3.*(-1.35e2./4.0)-cos(th5).*(4.05e2./4.0);t3.*(-1.35e2./4.0)];
