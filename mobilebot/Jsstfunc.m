function Jsst = Jsstfunc(in1)
%JSSTFUNC
%    JSST = JSSTFUNC(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    10-May-2018 12:26:44

th1 = in1(1,:);
th2 = in1(2,:);
th3 = in1(3,:);
th4 = in1(4,:);
th5 = in1(5,:);
t2 = sin(th4);
t3 = th3+3.0./1.0e1;
t4 = cos(th4);
t5 = th3.*1.0e1;
t6 = sin(th5);
t7 = t5-t6.*5.0+3.0;
t8 = t4.*th2;
Jsst = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,th2,-th1,0.0,0.0,0.0,1.0,t2.*t3,-t3.*t4,t8-t2.*th1,-t4,-t2,0.0,t2.*t7.*(1.0./1.0e1),t4.*t7.*(-1.0./1.0e1),t8+cos(th5).*(1.0./2.0)-t2.*th1,-t4,-t2,0.0],[6,6]);
