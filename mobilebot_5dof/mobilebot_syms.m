clear all
clc

syms th1 th2 th3 th4 th5 th6 real
syms th1d th2d th3d th4d th5d th6d real

th  = [th1  th2  th3 th4 th5 ]';
thd = [th1d th2d th3d th4d th5d ]';

n = 5;                    	% DOF

alumdens=2700;              % kg/m^3 for aluminum

g = 10;                    	% gravity, m/s^2 or N/kg

w = 1/10;                   % normal robot link width, meters

L0=0.3;                   	% link lengths (meters)
L1=0.5;
L2=0.5;

r0=L0/2;                    % distance to center of mass in meters
r1=L1/2;
r2=L2/2;

% Link home coordinates
gs10 = [eye(3) [0       0  0]';0 0 0 1];
gs20 = [eye(3) [0       0  0]';0 0 0 1];
gs30 = [eye(3) [-w/2    0 r0]';0 0 0 1];
gs40 = [eye(3) [w/2    r1 L0]';0 0 0 1];
gs50 = [eye(3) [w/2 L1+r2 L0]';0 0 0 1];

% Rectangular prisms ([xmin xmax ymin ymax zmin zmax]) in meters
RP = [  0   0   0   0   0   0;
        0   0   0   0   0   0;
      -2*w 2*w  -2*w 2*w -L0/2  L0/2;
      -w/2 w/2 -L1/2 L1/2  -w/2   w/2;
      -w/2 w/2 -L2/2 L2/2  -w/2   w/2]

% save dimensions etc. to a mat file to use in animation file
save('dimensions_4_3.mat','n','L0','L1','L2','r0','r1','r2','w','RP') 

% Calculate mass
m  = alumdens*(RP(:,2)-RP(:,1)).*(RP(:,4)-RP(:,3)).*(RP(:,6)-RP(:,5))

% Calculate moments of inertia
Ix = (m.*((RP(:,4)-RP(:,3)).^2+(RP(:,6)-RP(:,5)).^2)/12)  % using link dimensions
Iy = (m.*((RP(:,2)-RP(:,1)).^2+(RP(:,6)-RP(:,5)).^2)/12)  % using link dimensions
Iz = (m.*((RP(:,4)-RP(:,3)).^2+(RP(:,2)-RP(:,1)).^2)/12)  % using link dimensions

Ix = round(Ix*100)/100% rounding to 2 decimal place
Iy = round(Iy*100)/100 % rounding to 2 decimal place
Iz = round(Iz*100)/100 % rounding to 2 decimal place

% Link axis in home position
w1=[ 1 0 0]';
w2=[ 0 1 0]';
w3=[ 0 0 1]';
w4=[-1 0 0]';
w5=[-1 0 0]';

% Link axis locations in home positions
q3=[-w/2 0  0]';
q4=[0    0 L0]';
q5=[0   L1 L0]';

% Link inertia matrices
M1 = [m(1)*eye(3) zeros(3);zeros(3) [Ix(1) 0 0;0 Iy(1) 0;0 0 Iz(1)]]
M2 = [m(2)*eye(3) zeros(3);zeros(3) [Ix(2) 0 0;0 Iy(2) 0;0 0 Iz(2)]]
M3 = [m(3)*eye(3) zeros(3);zeros(3) [Ix(3) 0 0;0 Iy(3) 0;0 0 Iz(3)]]
M4 = [m(4)*eye(3) zeros(3);zeros(3) [Ix(4) 0 0;0 Iy(4) 0;0 0 Iz(4)]]
M5 = [m(5)*eye(3) zeros(3);zeros(3) [Ix(5) 0 0;0 Iy(5) 0;0 0 Iz(5)]]


% Twists
z1 = twistp(w1);
z2 = twistp(w2);
z3 = twistr(w3,q3);
z4 = twistr(w4,q4);
z5 = twistr(w5,q5);
z  = [z1 z2 z3 z4 z5];              % put in list for access later

z2_p = simplify(Ad(expt(z1,th1))*z2)
z3_p = simplify(Ad(expt(z1,th1)*expt(z2,th2))*z3)
z4_p = simplify(Ad(expt(z1,th1)*expt(z2,th2)*expt(z3,th3))*z4)
z5_p = simplify(Ad(expt(z1,th1)*expt(z2,th2)*expt(z3,th3)*expt(z4,th4))*z5)

Jsst = [z1 z2_p z3_p z4_p z5_p]

% Homogeneous transformation to each link frame
gs1th = simplify(expt(z1,th1)*gs10);
gs2th = simplify(expt(z1,th1)*expt(z2,th2)*gs20);
gs3th = simplify(expt(z1,th1)*expt(z2,th2)*expt(z3,th3)*gs30);
gs4th = simplify(expt(z1,th1)*expt(z2,th2)*expt(z3,th3)*expt(z4,th4)*gs40);
gs5th = simplify(expt(z1,th1)*expt(z2,th2)*expt(z3,th3)*expt(z4,th4)*expt(z5,th5)*gs50);



% Create the A matrices, put into matlag structure variable for later
for i = 1:n
    for j = 1:n
        if i<j
            A(i,j).M = zeros(6);
        elseif i==j
            A(i,j).M = eye(6);
        elseif i>j
            T=eye(4);
            for k = j+1:i
                T=T*expt(z(:,k),th(k));
            end
            A(i,j).M = Ad(ginv(T));
        end
    end
end

% Calculate the transformed link inertia matrices, put into matlab struct
% variable for later access
Mp(1).M = sym(Ad(ginv(gs10))'*M1*Ad(ginv(gs10)));
Mp(2).M = sym(Ad(ginv(gs20))'*M2*Ad(ginv(gs20)));
Mp(3).M = sym(Ad(ginv(gs30))'*M3*Ad(ginv(gs30)));
Mp(4).M = sym(Ad(ginv(gs40))'*M4*Ad(ginv(gs40)));
Mp(5).M = sym(Ad(ginv(gs50))'*M5*Ad(ginv(gs50)));


% Calculate the Inertia Matrix
M=sym(zeros(n,n));
for i = 1:n
    for j = 1:n
        M(i,j) = 0;
        for l = max(i,j):n
            M(i,j) = M(i,j)+z(:,i)'*A(l,i).M'*Mp(l).M*A(l,j).M*z(:,j);
        end
    end
end
M=simplify(M)

% Calculate the Coriolis Matrix
C=sym(zeros(n,n));
for i = 1:n
    for j = 1:n
        C(i,j)=0;
        for k = 1:n
            C(i,j)=C(i,j)+(Csub(i,j,k,n,z,A,Mp)+Csub(i,k,j,n,z,A,Mp)-Csub(k,j,i,n,z,A,Mp))*thd(k);
        end
        C(i,j)=C(i,j)/2;
    end
end
C=simplify(C)

% Calculate the external forces (gravity only, in this case)
cm1 = gs1th*[0 0 0 1]';  %center of mass of link 1 w.r.t spatial frame
cm2 = gs2th*[0 0 0 1]';
cm3 = gs3th*[0 0 0 1]';
cm4 = gs4th*[0 0 0 1]';
cm5 = gs5th*[0 0 0 1]';

V = g*[m(1) m(2) m(3) m(4) m(5)]*[cm1(3) cm2(3) cm3(3) cm4(3) cm5(3)]'   %summation of mgh for all 3 links
N = simplify(jacobian(V,th)')


%% STATE 
x = [th;thd]

%% create Matlab functions

matlabFunction(M,'file','Mfunc','vars',{x})
matlabFunction(C,'file','Cfunc','vars',{x})
matlabFunction(N,'file','Nfunc','vars',{x})

matlabFunction(gs1th,'file','gs1func','vars',{x})
matlabFunction(gs2th,'file','gs2func','vars',{x})
matlabFunction(gs3th,'file','gs3func','vars',{x})
matlabFunction(gs4th,'file','gs4func','vars',{x})
matlabFunction(gs5th,'file','gs5func','vars',{x})


matlabFunction(Jsst,'file','Jsstfunc','vars',{x})

%% create blocks directly for Simulink!

% new_system('Example4_3_Model_matlabFunctionBlock'); % only use if file is not open!
% open_system('Example4_3_Model_matlabFunctionBlock'); % only use if file is not open!

matlabFunctionBlock('mobilebot_Model/Mfunc',M,'vars',{x})
matlabFunctionBlock('mobilebot_Model/Cfunc',C,'vars',{x})
matlabFunctionBlock('mobilebot_Model/Nfunc',N,'vars',{x})

%matlabFunctionBlock('Example4_3_Model/Jsstfunc',Jsst,'vars',{x})

