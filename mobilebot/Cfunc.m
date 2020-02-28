function C = Cfunc(in1)
%CFUNC
%    C = CFUNC(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    10-May-2018 12:26:43

th4 = in1(4,:);
th5 = in1(5,:);
th6 = in1(6,:);
th4d = in1(10,:);
th5d = in1(11,:);
th6d = in1(12,:);
t2 = cos(th4);
t3 = cos(th5);
t4 = sin(th4);
t5 = cos(th6);
t6 = sin(th5);
t7 = sin(th6);
t8 = t2.*t3.*t7.*th4d.*(2.7e1./8.0);
t9 = t2.*t5.*t6.*th4d.*(2.7e1./8.0);
t10 = t3.*t4.*t5.*th5d.*(2.7e1./8.0);
t11 = t3.*t4.*t5.*th6d.*(2.7e1./8.0);
t12 = t3.*t4.*t7.*th4d.*(2.7e1./8.0);
t13 = t4.*t5.*t6.*th4d.*(2.7e1./8.0);
t14 = t2.*t6.*t7.*th5d.*(2.7e1./8.0);
t15 = t2.*t6.*t7.*th6d.*(2.7e1./8.0);
t16 = th5+th6;
t17 = sin(t16);
t18 = th5.*2.0;
t19 = th6.*2.0;
t20 = t18+t19;
t21 = sin(t20);
t22 = t18+th6;
t23 = sin(t22);
t24 = sin(t18);
t25 = cos(t16);
t26 = th5d+th6d;
t27 = t21.*3.3e1;
C = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2.*th4d.*(5.13e2./1.0e2)+t3.*t4.*th4d.*(8.1e1./8.0)+t2.*t6.*th5d.*(8.1e1./8.0)+t3.*t4.*t5.*th4d.*(2.7e1./8.0)+t2.*t3.*t7.*th5d.*(2.7e1./8.0)+t2.*t3.*t7.*th6d.*(2.7e1./8.0)+t2.*t5.*t6.*th5d.*(2.7e1./8.0)+t2.*t5.*t6.*th6d.*(2.7e1./8.0)-t4.*t6.*t7.*th4d.*(2.7e1./8.0),t4.*th4d.*(5.13e2./1.0e2)-t2.*t3.*th4d.*(8.1e1./8.0)+t4.*t6.*th5d.*(8.1e1./8.0)-t2.*t3.*t5.*th4d.*(2.7e1./8.0)+t2.*t6.*t7.*th4d.*(2.7e1./8.0)+t3.*t4.*t7.*th5d.*(2.7e1./8.0)+t3.*t4.*t7.*th6d.*(2.7e1./8.0)+t4.*t5.*t6.*th5d.*(2.7e1./8.0)+t4.*t5.*t6.*th6d.*(2.7e1./8.0),0.0,t7.*th6d.*(-2.7e1./3.2e1)-t21.*th5d.*5.56875e-1-t21.*th6d.*5.56875e-1-t23.*th5d.*(2.7e1./1.6e1)-t23.*th6d.*(2.7e1./3.2e1)-t24.*th5d.*2.244375,th4d.*(t23.*1.0e2+t24.*1.33e2+t27).*1.6875e-2,th4d.*(t7.*5.0e1+t23.*5.0e1+t27).*1.6875e-2,t8+t9+t10+t11+t2.*t6.*th4d.*(8.1e1./8.0)+t3.*t4.*th5d.*(8.1e1./8.0)-t4.*t6.*t7.*th5d.*(2.7e1./8.0)-t4.*t6.*t7.*th6d.*(2.7e1./8.0),t12+t13+t14+t15-t2.*t3.*th5d.*(8.1e1./8.0)+t4.*t6.*th4d.*(8.1e1./8.0)-t2.*t3.*t5.*th5d.*(2.7e1./8.0)-t2.*t3.*t5.*th6d.*(2.7e1./8.0),t6.*th5d.*(8.1e1./8.0)+t17.*th5d.*(2.7e1./8.0)+t17.*th6d.*(2.7e1./8.0),t3.*th5d.*(-8.1e1./1.6e2)-t21.*th4d.*5.56875e-1-t23.*th4d.*(2.7e1./1.6e1)-t24.*th4d.*2.244375-t25.*th5d.*(2.7e1./1.6e2)-t25.*th6d.*(2.7e1./1.6e2),t7.*th6d.*(-2.7e1./1.6e1),t7.*th5d.*(2.7e1./1.6e1),t8+t9+t10+t11-t4.*t6.*t7.*th5d.*(2.7e1./8.0)-t4.*t6.*t7.*th6d.*(2.7e1./8.0),t12+t13+t14+t15-t2.*t3.*t5.*th5d.*(2.7e1./8.0)-t2.*t3.*t5.*th6d.*(2.7e1./8.0),t17.*t26.*(2.7e1./8.0),t7.*th4d.*(-2.7e1./3.2e1)-t21.*th4d.*5.56875e-1-t23.*th4d.*(2.7e1./3.2e1)-t25.*th5d.*(2.7e1./1.6e2)-t25.*th6d.*(2.7e1./1.6e2),t7.*t26.*(-2.7e1./1.6e1),0.0],[6,6]);