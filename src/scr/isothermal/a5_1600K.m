
clear;clc;
time0=clock();
format long;
%%%%%%%%%%%%%%%%%%%%%%%%%%%paraments%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%fe-27.35wt%cr-17.15wt%ni
Nx = 300;
Ny = 300;
NxNy= Nx*Ny;
N=0.3333;

dx=1e-8;
dy=1e-8;
dt=1e-10;
ddt=dx^2/(6*dx);
tol=1e-12;


nstep =1e16;
nprint =10000;

Vm=7.1e-6;                 
sigma12=1;  
sigma13=0.49;  
sigma23=1; 
sigma123=15*sigma12;

 

epsilon=2*dx;          % m
ttime=0.0;
x=0.02;%0.01           
%%%%%%%%%%%%%%%%%%%%%%%%%diffusion coefficient%%%%%%%%%%%%%%%%%%%%%%%
D11 =  +2.08556E-11;%%m^2/s    TCFE12
D12 =  +1.45564E-11;
D13 = +7.84621E-12;
D14 =  +1.71842E-11;

D31 =  +2.77638E-09;
D32 =  +1.74294E-10;
D33 =  +1.02311E-10;
D34 =  +2.64079E-09;


% D11=0.1*D31;
% D12=0.1*D32;
% D13=0.1*D33;
% D14=0.1*D34;
% 
% D21=0.1*D31;
% D22=0.1*D32;
% D23=0.1*D33;
% D24=0.1*D34;

D21=9.06642E-14;            %cc
D22= -4.81719E-14;            %cmn
D23=-1.33761E-14;            %MN C
D24=1.37852E-13;            %MN MN


G1=-97293.72;   %bcc J/mol      TCFE12
G2=-103808.87;
G3=-98451.158;   %liq

miui =-88221.545;             %CR J/mol       TCFE12
miuj = -114065.5;             %Ni
miuk = -97523.374;             %Fe

dG1Dki=-23682.707;
dG1Dkj=-25386.755;

dG3Dki=-29460.693;
dG3Dkj=-26599.986;

dG2Dki=0;
dG2Dkj=0;



    dG1Dii= 38850.4-dG1Dki;       %TCFE12      
    dG1Dij= 888.84973-dG1Dkj;
    dG1Dji= 2592.8978-dG1Dki;
    dG1Djj= 81509.414-dG1Dkj;

    dG3Dii= 54417.866-dG3Dki;
    dG3Dij= 3024.0409-dG3Dkj;
    dG3Dji= 163.33369-dG3Dki;
    dG3Djj= 56255.437-dG3Dkj;
   
    dG2Dii =50892.198;
    dG2Dij =-4165.8759;
    dG2Dji =-6276.186;
    dG2Djj =58120.011;
   
        c21=0.31719646;       %C
       c22=0.25497879;     %Mn
       c23=1-c21-c22;
 

c11=0.31106007;   
c12=0.16102975;    
c13=1-c11-c12;
        
c31=0.27222039;       
c32=0.22371163;     
c33=1-c31-c32;

miu=(miui-miuk)/Vm;
miu2=(miuj-miuk)/Vm;

ddt2=0.25*dx.^2/D13;
%%%%%%%%%%%%%%%%%%%%%%%%%%fitting coefficient%%%%%%%%%%%%%%%%%%%%%%%%%%%
C1=0.5*dG1Dji/Vm;                                                 %j/m3
A1=0.5*dG1Dii/Vm-C1;
B1=0.5*dG1Djj/Vm-C1;
O1 = miu-2*A1*c11+2*C1*c13;          
P1 = miu2-2*B1*c12+2*C1*c13;
Q1 = G1/Vm-(A1*c11.^2+B1*c12.^2+C1*c13.^2+O1*c11+P1*c12);
        
   
C2=0.5*dG2Dji/Vm;                                                 %j/m3
A2=0.5*dG2Dii/Vm-C2;
B2=0.5*dG2Djj/Vm-C2;
O2 = miu-2*A2*c21+2*C2*c23;            
P2 = miu2-2*B2*c22+2*C2*c23;
Q2 = G2/Vm-(A2*c21.^2+B2*c22.^2+C2.*c23.^2+O2.*c21+P2.*c22);
        

C3=0.5*dG3Dji/Vm;                                                 %j/m3
A3=0.5*dG3Dii/Vm-C3;
B3=0.5*dG3Djj/Vm-C3;
O3 = miu-2*A3*c31+2*C3*c33;          
P3 = miu2-2*B3*c32+2*C3*c33;
Q3 = G3/Vm-(A3*c31.^2+B3*c32.^2+C3*c33.^2+O3*c31+P3*c32);
        

dCom1_dmu1=0.5/C1/((A1+C1)/C1-C1/(B1+C1));
dCom1_dmu2=-1/2/(B1+C1)/((A1+C1)/C1-C1/(B1+C1));
dCom1_dmu3=1/2/(A1+C1)/(C1/(A1+C1)-(B1+C1)/C1);
dCom1_dmu4=-1/2/C1/(C1/(A1+C1)-(B1+C1)/C1);
  
         
dCom2_dmu1=0.5/C2/((A2+C2)/C2-C2/(B2+C2));
dCom2_dmu2=-1/2/(B2+C2)/((A2+C2)/C2-C2/(B2+C2));
dCom2_dmu3=1/2/(A2+C2)/(C2/(A2+C2)-(B2+C2)/C2);
dCom2_dmu4=-1/2/C2/(C2/(A2+C2)-(B2+C2)/C2);
                
        
dCom3_dmu1=0.5/C3/((A3+C3)/C3-C3/(B3+C3));
dCom3_dmu2=-1/2/(B3+C3)/((A3+C3)/C3-C3/(B3+C3));
dCom3_dmu3=1/2/(A3+C3)/(C3/(A3+C3)-(B3+C3)/C3);
dCom3_dmu4=-1/2/C3/(C3/(A3+C3)-(B3+C3)/C3);
% tau13=7.72e9;
% tau23=7.1e9;
% tau12=7.1e15;   
% tau13 = epsilon * 0.182 * (((c11 - c31).^2 / (D31 * dCom3_dmu1))+(2*(c31 - c11)*(c32 - c12)/(D32 * dCom3_dmu2))+((c32 - c12).^2/(D34 * dCom3_dmu4)));
% tau131=((c11 - c31).^2 / (D31 * dCom3_dmu1));
% tau132=(2*(c31 - c11)*(c32 - c12)/(D32 * dCom3_dmu2));
% tau133=((c32 - c12).^2/(D34 * dCom3_dmu4));
% tau23 = epsilon * 0.182 * (((c21 - c31).^2 / (D31 * dCom3_dmu1))+(2*(c31 - c21)*(c32 - c22)/(D32 * dCom3_dmu2))+((c32 - c22).^2/(D34 * dCom3_dmu4)));
% tau12 = epsilon * 0.182 * (((c11 - c21).^2 / (D11 * dCom1_dmu1))+(2*(c11 - c21)*(c21 - c11)/(D11 * dCom1_dmu3))+((c12 - c22).^2/(D12 * dCom1_dmu4)));
% 
% tau13 = epsilon * 0.223 * (c11 - c31)^2 / (D31 * dCom3_dmu1);
% tau23 = epsilon * 0.223 * (c21 - c31)^2 / (D31 * dCom3_dmu1);
% tau12= epsilon * 0.223 * (c21 - c11)^2 / (D21 * dCom2_dmu1);
% tau12=tau23;
Tau123_1=5.44e7;
Tau123_2=1.12e7;
Tau123_3=5.44e7;
% Tau123=(Tau123_1+Tau123_2+Tau123_3)/3;
Tau123=4e10;
A_D1=[D11*dCom1_dmu1,D12*dCom1_dmu2;D13*dCom1_dmu3,D14*dCom1_dmu4];
dA_D1=inv(A_D1);
B_D1=[D11*dCom1_dmu1,D12*dCom1_dmu2;D13*dCom1_dmu3,D14*dCom1_dmu4];
dB_D1=inv(B_D1);
tau13 = epsilon * 0.223 * [c31-c11,c32-c12]*inv(A_D1)*[c31-c11;c32-c12];
tau23 = epsilon * 0.223 * [c31-c21,c32-c22]*inv(A_D1)*[c31-c21;c32-c22];     
tau12 = epsilon * 0.223 * [c21-c11,c22-c12]*inv(B_D1)*[c21-c11;c22-c12]; 
%tau12=tau23;
tau13=3.6e9;


Phi1=zeros(Nx,Ny);
Phi2=zeros(Nx,Ny);
Phi3=zeros(Nx,Ny);
Phi1_dot = zeros(Nx, Ny);
Phi2_dot = zeros(Nx, Ny);
Phi3_dot = zeros(Nx, Ny);
Phi_dot= zeros(Nx, Ny);

gr_po1=zeros(Nx,Ny);
gr_po3=zeros(Nx,Ny);
gr_po2=zeros(Nx,Ny);


Com11=zeros(Nx,Ny);
Com13=zeros(Nx,Ny);
Com12=zeros(Nx,Ny);
Com21=zeros(Nx,Ny);
Com22=zeros(Nx,Ny);
Com23=zeros(Nx,Ny);
Com31=zeros(Nx,Ny);
Com32=zeros(Nx,Ny);
Com33=zeros(Nx,Ny);

lap_Phi1=zeros(Nx,Ny);
lap_Phi2=zeros(Nx,Ny);
lap_Phi3=zeros(Nx,Ny);

grad_x_Phi1=zeros(Nx,Ny);
grad_y_Phi1=zeros(Nx,Ny);
grad_x_Phi2=zeros(Nx,Ny);
grad_y_Phi2=zeros(Nx,Ny);
grad_x_Phi3=zeros(Nx,Ny);
grad_y_Phi3=zeros(Nx,Ny);

g_Phi11=zeros(Nx,Ny);
g_Phi22=zeros(Nx,Ny);
g_Phi33=zeros(Nx,Ny);
g_Phi12=zeros(Nx,Ny);
g_Phi13=zeros(Nx,Ny);
g_Phi23=zeros(Nx,Ny);

grad_J11=zeros(Nx,Ny);
grad_J12=zeros(Nx,Ny);
grad_J21=zeros(Nx,Ny);
grad_J22=zeros(Nx,Ny);

GG1=zeros(Nx,Ny);
GG2=zeros(Nx,Ny);
GG3=zeros(Nx,Ny);

termg33_2=zeros(Nx,Ny);
termh33_2=zeros(Nx,Ny);
termi33_2=zeros(Nx,Ny);

W12=zeros(Nx,Ny);
W21=zeros(Nx,Ny);
W13=zeros(Nx,Ny);
W31=zeros(Nx,Ny);
W23=zeros(Nx,Ny);
W32=zeros(Nx,Ny);


hp_Phi1=zeros(Nx,Ny);
hp_Phi2=zeros(Nx,Ny);
hp_Phi3=zeros(Nx,Ny);
dh_Phi1=zeros(Nx,Ny);
dh_Phi2=zeros(Nx,Ny);
dh_Phi3=zeros(Nx,Ny);

J11_1=zeros(Nx,Ny);
J21_1=zeros(Nx,Ny);
J12_1=zeros(Nx,Ny);
J22_1=zeros(Nx,Ny);

x2_J11=zeros(Nx,Ny);
x1_J11=zeros(Nx,Ny);
y2_J11=zeros(Nx,Ny);
y1_J11=zeros(Nx,Ny);
x2_J21=zeros(Nx,Ny);
x1_J21=zeros(Nx,Ny);
y2_J21=zeros(Nx,Ny);
y1_J21=zeros(Nx,Ny);

x2_J11_1=zeros(Nx,Ny);
x1_J11_1=zeros(Nx,Ny);
y2_J11_1=zeros(Nx,Ny);
y1_J11_1=zeros(Nx,Ny);
x2_J21_1=zeros(Nx,Ny);
x1_J21_1=zeros(Nx,Ny);
y2_J21_1=zeros(Nx,Ny);
y1_J21_1=zeros(Nx,Ny);

x2_J12=zeros(Nx,Ny);
x1_J12=zeros(Nx,Ny);
y2_J12=zeros(Nx,Ny);
y1_J12=zeros(Nx,Ny);
x2_J22=zeros(Nx,Ny);
x1_J22=zeros(Nx,Ny);
y2_J22=zeros(Nx,Ny);
y1_J22=zeros(Nx,Ny);

x2_J12_1=zeros(Nx,Ny);
x1_J12_1=zeros(Nx,Ny);
y2_J12_1=zeros(Nx,Ny);
y1_J12_1=zeros(Nx,Ny);
x2_J22_1=zeros(Nx,Ny);
x1_J22_1=zeros(Nx,Ny);
y2_J22_1=zeros(Nx,Ny);
y1_J22_1=zeros(Nx,Ny);



numer_x1_J11=zeros(Nx,Ny);
numer_x2_J11=zeros(Nx,Ny);
numer_y1_J11=zeros(Nx,Ny);
numer_y2_J11=zeros(Nx,Ny);
numer_x1_J21=zeros(Nx,Ny);
numer_x2_J21=zeros(Nx,Ny);
numer_y1_J21=zeros(Nx,Ny);
numer_y2_J21=zeros(Nx,Ny);
numer_x1_J12=zeros(Nx,Ny);
numer_x2_J12=zeros(Nx,Ny);
numer_y1_J12=zeros(Nx,Ny);
numer_y2_J12=zeros(Nx,Ny);
numer_x1_J22=zeros(Nx,Ny);
numer_x2_J22=zeros(Nx,Ny);
numer_y1_J22=zeros(Nx,Ny);
numer_y2_J22=zeros(Nx,Ny);

denom_x1_J11=zeros(Nx,Ny);
denom_x2_J11=zeros(Nx,Ny);
denom_y1_J11=zeros(Nx,Ny);
denom_y2_J11=zeros(Nx,Ny);
denom_x1_J21=zeros(Nx,Ny);
denom_x2_J21=zeros(Nx,Ny);
denom_y1_J21=zeros(Nx,Ny);
denom_y2_J21=zeros(Nx,Ny);
denom_x1_J12=zeros(Nx,Ny);
denom_x2_J12=zeros(Nx,Ny);
denom_y1_J12=zeros(Nx,Ny);
denom_y2_J12=zeros(Nx,Ny);
denom_x1_J22=zeros(Nx,Ny);
denom_x2_J22=zeros(Nx,Ny);
denom_y1_J22=zeros(Nx,Ny);
denom_y2_J22=zeros(Nx,Ny);

denom_x1_J11w_x=zeros(Nx,Ny);
denom_x1_J11w_y=zeros(Nx,Ny);
denom_x1_J21w_x=zeros(Nx,Ny);
denom_x1_J21w_y=zeros(Nx,Ny);
denom_x2_J11e_x=zeros(Nx,Ny);
denom_x2_J11e_y=zeros(Nx,Ny);
denom_x2_J21e_x=zeros(Nx,Ny);
denom_x2_J21e_y=zeros(Nx,Ny);
denom_y1_J11s_x=zeros(Nx,Ny);
denom_y1_J11s_y=zeros(Nx,Ny);
denom_y2_J21s_x=zeros(Nx,Ny);
denom_y2_J21s_y=zeros(Nx,Ny);
denom_y2_J11n_x=zeros(Nx,Ny);
denom_y2_J11n_y=zeros(Nx,Ny);
denom_y2_J21n_x=zeros(Nx,Ny);
denom_y2_J21n_y=zeros(Nx,Ny);
denom_y1_J21s_x=zeros(Nx,Ny);
denom_y1_J21s_y=zeros(Nx,Ny);
denom_x2_J12e_x=zeros(Nx,Ny);
denom_x2_J12e_y=zeros(Nx,Ny);
denom_x1_J12w_x=zeros(Nx,Ny);
denom_x1_J12w_y=zeros(Nx,Ny);
denom_y2_J12n_x=zeros(Nx,Ny);
denom_y2_J12n_y=zeros(Nx,Ny);
denom_y1_J12s_x=zeros(Nx,Ny);
denom_y1_J12s_y=zeros(Nx,Ny);
denom_x2_J22e_x=zeros(Nx,Ny);
denom_x2_J22e_y=zeros(Nx,Ny);
denom_x1_J22w_x=zeros(Nx,Ny);
denom_x1_J22w_y=zeros(Nx,Ny);
denom_y2_J22n_x=zeros(Nx,Ny);
denom_y2_J22n_y=zeros(Nx,Ny);
denom_y1_J22s_x=zeros(Nx,Ny);
denom_y1_J22s_y=zeros(Nx,Ny);









J11=zeros(Nx,Ny);
J12=zeros(Nx,Ny);
J21=zeros(Nx,Ny);
J22=zeros(Nx,Ny);

termg=zeros(Nx,Ny);
termh=zeros(Nx,Ny);
termi=zeros(Nx,Ny);

 Term3=zeros(Nx,Ny);
 Term2=zeros(Nx,Ny);
 Term1=zeros(Nx,Ny);

a=zeros(Nx,Ny);
b=zeros(Nx,Ny);
c=zeros(Nx,Ny);
d=zeros(Nx,Ny);

fenzi1=zeros(Nx,Ny);
fenzi2=zeros(Nx,Ny);
fenmu=zeros(Nx,Ny);
fenzi12=zeros(Nx,Ny);
fenzi22=zeros(Nx,Ny);
miuu1=zeros(Nx,Ny);
miuu2=zeros(Nx,Ny);
m=zeros(Nx,Ny);
m1=zeros(Nx,Ny);
m2=zeros(Nx,Ny);
m3=zeros(Nx,Ny);
m4=zeros(Nx,Ny);

%%%%%%%%%%%%%%%%%% Phase initial%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:Nx
   for j = 1:Ny
      xx0 =sqrt((i-75)^2 +(j-100)^2);
        xx1 =sqrt((i-225)^2 +(j-50)^2);
        xx2 =sqrt((i-150)^2 +(j-75)^2);
        xx3 =sqrt((i-120)^2 +(j-150)^2);
        xx4 =sqrt((i-200)^2 +(j-125)^2);
        xx5 =sqrt((i-50)^2 +(j-150)^2);
        xx6 =sqrt((i-225)^2 +(j-175)^2);
        xx7 =sqrt((i-150)^2 +(j-200)^2);
        xx8 =sqrt((i-50)^2 +(j-225)^2);
        xx9 =sqrt((i-170)^2 +(j-265)^2);
        


               if(xx0<12||xx1<12||xx2<10||xx3<18||xx4<18||xx5<18||xx6<15||xx7<15||xx8<19||xx9<20)
                 Phi1(i,j)=1.0;
                 Phi2(i,j)=0.0;  
                 Phi3(i,j)=0.0;
              % elseif(xx1<10||xx2<5||xx3<7)
              %    Phi1(i,j)=0.0;
              %    Phi2(i,j)=1.0;  
              %    Phi3(i,j)=0.0;    
              else
                 Phi1(i,j)=0.0;
                 Phi2(i,j)=0.0;  
                 Phi3(i,j)=1.0;  
              end

              

            
              
         
    end
end

miu(1:Nx,1:Ny)=(miui-miuk)/Vm;     
miu2(1:Nx,1:Ny)=(miuj-miuk)/Vm;    
tic;

for istep =1:nstep
   ttime=ttime+dt;

        Com11(1:Nx,1:Ny)=((miu(1:Nx,1:Ny)-O1+2*C1)/2/C1-(miu2(1:Nx,1:Ny)-P1+2*C1)/(B1+C1)/2)/((A1+C1)/C1-C1/(B1+C1));
        Com21(1:Nx,1:Ny)=((miu(1:Nx,1:Ny)-O2+2*C2)/2/C2-(miu2(1:Nx,1:Ny)-P2+2*C2)/(B2+C2)/2)/((A2+C2)/C2-C2/(B2+C2));
        Com31(1:Nx,1:Ny)=((miu(1:Nx,1:Ny)-O3+2*C3)/2/C3-(miu2(1:Nx,1:Ny)-P3+2*C3)/(B3+C3)/2)/((A3+C3)/C3-C3/(B3+C3));
        Com12(1:Nx,1:Ny)=((miu(1:Nx,1:Ny)-O1+2*C1)/2/(A1+C1)-(miu2(1:Nx,1:Ny)-P1+2*C1)/2/C1)/(C1/(A1+C1)-(B1+C1)/C1);
        Com22(1:Nx,1:Ny)=((miu(1:Nx,1:Ny)-O2+2*C2)/2/(A2+C2)-(miu2(1:Nx,1:Ny)-P2+2*C2)/2/C2)/(C2/(A2+C2)-(B2+C2)/C2);
        Com32(1:Nx,1:Ny)=((miu(1:Nx,1:Ny)-O3+2*C3)/2/(A3+C3)-(miu2(1:Nx,1:Ny)-P3+2*C3)/2/C3)/(C3/(A3+C3)-(B3+C3)/C3);
        Com13(1:Nx,1:Ny)=(1-Com11(1:Nx,1:Ny)-Com12(1:Nx,1:Ny));
        Com23(1:Nx,1:Ny)=(1-Com21(1:Nx,1:Ny)-Com22(1:Nx,1:Ny));
        Com33(1:Nx,1:Ny)=(1-Com31(1:Nx,1:Ny)-Com32(1:Nx,1:Ny));

        % Com3=1-Com1(1:Nx,1:Ny)-Com2(1:Nx,1:Ny);
        
        
 %%%%%%%%%%%%%%%Update of Phase Field
        gr_po1(1:Nx,1:Ny) = Q1-A1.*Com11(1:Nx,1:Ny).^2-B1.*Com12(1:Nx,1:Ny).^2+C1*Com13(1:Nx,1:Ny).^2+...
                            2*C1*Com11(1:Nx,1:Ny).*Com13(1:Nx,1:Ny)+2*C1*Com13(1:Nx,1:Ny).*Com12(1:Nx,1:Ny);
        gr_po2(1:Nx,1:Ny) = Q2-A2.*Com21(1:Nx,1:Ny).^2-B2.*Com22(1:Nx,1:Ny).^2+C2*Com23(1:Nx,1:Ny).^2+...
                            2*C2*Com21(1:Nx,1:Ny).*Com23(1:Nx,1:Ny)+2*C2*Com23(1:Nx,1:Ny).*Com22(1:Nx,1:Ny);
        gr_po3(1:Nx,1:Ny) = Q3-A3.*Com31(1:Nx,1:Ny).^2-B3.*Com32(1:Nx,1:Ny).^2+C3*Com33(1:Nx,1:Ny).^2+...
                            2*C3*Com31(1:Nx,1:Ny).*Com33(1:Nx,1:Ny)+2*C3*Com33(1:Nx,1:Ny).*Com32(1:Nx,1:Ny);

        g1_3=gr_po1(1:Nx,1:Ny)-gr_po3(1:Nx,1:Ny);

       

for i=1:Nx 
    for j=1:Ny
         jp=j+1;
         jm=j-1;
         ip=i+1;
         im=i-1;
        if(im == 0)
           im=Nx;
        end
        if(ip == (Nx+1))
           ip=1;
        end
       if(jm == 0) 
           jm = 1;
        end
        if(jp == (Ny+1))
           jp=Ny;
        end
        % if(jm == 0) 
        %    jm = Ny;
        % end
        % if(jp == (Ny+1))
        %    jp=1;
        % end



        lap_Phi1(i,j)=(Phi1(ip,j)+Phi1(im,j)+Phi1(i,jm)+Phi1(i,jp) ...
                     +0.5*(Phi1(ip,jm)+Phi1(im,jm)+Phi1(ip,jp)+Phi1(im,jp)) ...
                     -6*Phi1(i,j))/((dx^2));
        lap_Phi2(i,j)=(Phi2(ip,j)+Phi2(im,j)+Phi2(i,jm)+Phi2(i,jp) ...
                     +0.5*(Phi2(ip,jm)+Phi2(im,jm)+Phi2(ip,jp)+Phi2(im,jp)) ...
                     -6*Phi2(i,j))/((dx^2));
        lap_Phi3(i,j)=(Phi3(ip,j)+Phi3(im,j)+Phi3(i,jm)+Phi3(i,jp) ...
                     +0.5*(Phi3(ip,jm)+Phi3(im,jm)+Phi3(ip,jp)+Phi3(im,jp)) ...
                     -6*Phi3(i,j))/((dx^2));


        grad_x_Phi1(i,j)=(Phi1(ip,j)-Phi1(im,j))/(2*dx);
        grad_y_Phi1(i,j)=(Phi1(i,jp)-Phi1(i,jm))/(2*dy);
        grad_x_Phi2(i,j)=(Phi2(ip,j)-Phi2(im,j))/(2*dx);
        grad_y_Phi2(i,j)=(Phi2(i,jp)-Phi2(i,jm))/(2*dy);
        grad_x_Phi3(i,j)=(Phi3(ip,j)-Phi3(im,j))/(2*dx);
        grad_y_Phi3(i,j)=(Phi3(i,jp)-Phi3(i,jm))/(2*dy);

   end
  end
     
    g_Phi11(1:Nx,1:Ny)=(grad_x_Phi1(1:Nx,1:Ny).^2)+(grad_y_Phi1(1:Nx,1:Ny).^2);
    g_Phi22(1:Nx,1:Ny)=(grad_x_Phi2(1:Nx,1:Ny).^2)+(grad_y_Phi2(1:Nx,1:Ny).^2);
    g_Phi33(1:Nx,1:Ny)=(grad_x_Phi3(1:Nx,1:Ny).^2)+(grad_y_Phi3(1:Nx,1:Ny).^2);

    g_Phi12(1:Nx,1:Ny)=grad_x_Phi1(1:Nx,1:Ny).*grad_x_Phi2(1:Nx,1:Ny)+grad_y_Phi1(1:Nx,1:Ny).*grad_y_Phi2(1:Nx,1:Ny);
    g_Phi13(1:Nx,1:Ny)=grad_x_Phi1(1:Nx,1:Ny).*grad_x_Phi3(1:Nx,1:Ny)+grad_y_Phi1(1:Nx,1:Ny).*grad_y_Phi3(1:Nx,1:Ny);
    g_Phi23(1:Nx,1:Ny)=grad_x_Phi2(1:Nx,1:Ny).*grad_x_Phi3(1:Nx,1:Ny)+grad_y_Phi2(1:Nx,1:Ny).*grad_y_Phi3(1:Nx,1:Ny);

    % dh_Phi1(1:Nx,1:Ny)=6*Phi1(1:Nx,1:Ny).*(1-Phi1(1:Nx,1:Ny));
    for i = 1 : Nx
    for j = 1 : Ny

      
     
    if((Phi1(i, j) > 0 && Phi1(i, j) < 1) && (Phi2(i, j) > 0 && Phi2(i, j) < 1) && (Phi3(i, j) > 0 && Phi3(i, j) < 1))
         dh_Phi1(i,j)=6*Phi1(i,j)*(1-Phi1(i,j))+2*Phi2(i,j)*Phi3(i,j);
         dh_Phi2(i,j)=6*Phi2(i,j)*(1-Phi2(i,j))+2*Phi1(i,j)*Phi3(i,j); 
         dh_Phi3(i,j)=6*Phi3(i,j)*(1-Phi3(i,j))+2*Phi1(i,j)*Phi2(i,j);

       

    else
    dh_Phi1(i,j)=6*Phi1(i,j)*(1-Phi1(i,j));
    dh_Phi2(i,j)=6*Phi2(i,j)*(1-Phi2(i,j)); 
    dh_Phi3(i,j)=6*Phi3(i,j)*(1-Phi3(i,j));

   
        
     end
    end
    end



    

    termg(1:Nx,1:Ny)=-epsilon*sigma12*(4*(Phi1(1:Nx,1:Ny).*g_Phi22(1:Nx,1:Ny)-Phi2(1:Nx,1:Ny).*g_Phi12(1:Nx,1:Ny))+2*Phi2(1:Nx,1:Ny).*(Phi1(1:Nx,1:Ny).*lap_Phi2(1:Nx,1:Ny)-Phi2(1:Nx,1:Ny).*lap_Phi1(1:Nx,1:Ny)))-(16/(pi^2))*(sigma12/epsilon)*(Phi2(1:Nx,1:Ny))-...
                      epsilon*sigma13*(4*(Phi1(1:Nx,1:Ny).*g_Phi33(1:Nx,1:Ny)-Phi3(1:Nx,1:Ny).*g_Phi13(1:Nx,1:Ny))+2*Phi3(1:Nx,1:Ny).*(Phi1(1:Nx,1:Ny).*lap_Phi3(1:Nx,1:Ny)-Phi3(1:Nx,1:Ny).*lap_Phi1(1:Nx,1:Ny)))-(16/(pi^2))*(sigma13/epsilon)*(Phi3(1:Nx,1:Ny))-(gr_po1(1:Nx,1:Ny).*dh_Phi1(1:Nx,1:Ny));
    termh(1:Nx,1:Ny)=-epsilon*sigma12*(4*(Phi2(1:Nx,1:Ny).*g_Phi11(1:Nx,1:Ny)-Phi1(1:Nx,1:Ny).*g_Phi12(1:Nx,1:Ny))+2*Phi1(1:Nx,1:Ny).*(Phi2(1:Nx,1:Ny).*lap_Phi1(1:Nx,1:Ny)-Phi1(1:Nx,1:Ny).*lap_Phi2(1:Nx,1:Ny)))-(16/(pi^2))*(sigma12/epsilon)*(Phi1(1:Nx,1:Ny))-...
                      epsilon*sigma23*(4*(Phi2(1:Nx,1:Ny).*g_Phi33(1:Nx,1:Ny)-Phi3(1:Nx,1:Ny).*g_Phi23(1:Nx,1:Ny))+2*Phi3(1:Nx,1:Ny).*(Phi2(1:Nx,1:Ny).*lap_Phi3(1:Nx,1:Ny)-Phi3(1:Nx,1:Ny).*lap_Phi2(1:Nx,1:Ny)))-(16/(pi^2))*(sigma23/epsilon)*(Phi3(1:Nx,1:Ny))-(gr_po2(1:Nx,1:Ny).*dh_Phi2(1:Nx,1:Ny));
    termi(1:Nx,1:Ny)=-epsilon*sigma13*(4*(Phi3(1:Nx,1:Ny).*g_Phi11(1:Nx,1:Ny)-Phi1(1:Nx,1:Ny).*g_Phi13(1:Nx,1:Ny))+2*Phi1(1:Nx,1:Ny).*(Phi3(1:Nx,1:Ny).*lap_Phi1(1:Nx,1:Ny)-Phi1(1:Nx,1:Ny).*lap_Phi3(1:Nx,1:Ny)))-(16/(pi^2))*(sigma13/epsilon)*(Phi1(1:Nx,1:Ny))-...
                      epsilon*sigma23*(4*(Phi3(1:Nx,1:Ny).*g_Phi22(1:Nx,1:Ny)-Phi2(1:Nx,1:Ny).*g_Phi23(1:Nx,1:Ny))+2*Phi2(1:Nx,1:Ny).*(Phi3(1:Nx,1:Ny).*lap_Phi2(1:Nx,1:Ny)-Phi2(1:Nx,1:Ny).*lap_Phi3(1:Nx,1:Ny)))-(16/(pi^2))*(sigma23/epsilon)*(Phi2(1:Nx,1:Ny))-(gr_po3(1:Nx,1:Ny).*dh_Phi3(1:Nx,1:Ny));
    %%%%2025.3.27 18:30  改双势阱
    termg33(1:Nx,1:Ny)=epsilon*(sigma13*(4*(g_Phi13(1:Nx,1:Ny).*Phi3(1:Nx,1:Ny)-Phi1(1:Nx,1:Ny).*g_Phi33(1:Nx,1:Ny))+2*Phi3(1:Nx,1:Ny).*(Phi3(1:Nx,1:Ny).*lap_Phi1(1:Nx,1:Ny)-Phi1(1:Nx,1:Ny).*lap_Phi3(1:Nx,1:Ny))) ...
                                +sigma12*(4*(g_Phi12(1:Nx,1:Ny).*Phi2(1:Nx,1:Ny)-Phi1(1:Nx,1:Ny).*g_Phi22(1:Nx,1:Ny))+2*Phi2(1:Nx,1:Ny).*(Phi2(1:Nx,1:Ny).*lap_Phi1(1:Nx,1:Ny)-Phi1(1:Nx,1:Ny).*lap_Phi2(1:Nx,1:Ny))));
                                

    termg33_2(1:Nx,1:Ny)=(16/(pi^2))*(sigma13*Phi3(1:Nx,1:Ny)+sigma12*Phi2(1:Nx,1:Ny)+sigma123*Phi2(1:Nx,1:Ny).*Phi3(1:Nx,1:Ny))/epsilon;

    termh33(1:Nx,1:Ny)=epsilon*(sigma23*(4*(g_Phi23(1:Nx,1:Ny).*Phi3(1:Nx,1:Ny)-Phi2(1:Nx,1:Ny).*g_Phi33(1:Nx,1:Ny))+2*Phi3(1:Nx,1:Ny).*(Phi3(1:Nx,1:Ny).*lap_Phi2(1:Nx,1:Ny)-Phi2(1:Nx,1:Ny).*lap_Phi3(1:Nx,1:Ny))) ...
                                +sigma12*(4*(g_Phi12(1:Nx,1:Ny).*Phi1(1:Nx,1:Ny)-Phi2(1:Nx,1:Ny).*g_Phi11(1:Nx,1:Ny))+2*Phi1(1:Nx,1:Ny).*(Phi1(1:Nx,1:Ny).*lap_Phi2(1:Nx,1:Ny)-Phi2(1:Nx,1:Ny).*lap_Phi1(1:Nx,1:Ny))));
                                

    termh33_2(1:Nx,1:Ny)=(16/(pi^2))*(sigma23*Phi3(1:Nx,1:Ny)+sigma12*Phi1(1:Nx,1:Ny)+sigma123*Phi1(1:Nx,1:Ny).*Phi3(1:Nx,1:Ny))/epsilon;
    
    termi33(1:Nx,1:Ny)=epsilon*(sigma13*(4*(g_Phi13(1:Nx,1:Ny).*Phi1(1:Nx,1:Ny)-Phi3(1:Nx,1:Ny).*g_Phi11(1:Nx,1:Ny))+2*Phi1(1:Nx,1:Ny).*(Phi1(1:Nx,1:Ny).*lap_Phi3(1:Nx,1:Ny)-Phi3(1:Nx,1:Ny).*lap_Phi1(1:Nx,1:Ny))) ...
                                +sigma23*(4*(g_Phi23(1:Nx,1:Ny).*Phi2(1:Nx,1:Ny)-Phi3(1:Nx,1:Ny).*g_Phi22(1:Nx,1:Ny))+2*Phi2(1:Nx,1:Ny).*(Phi2(1:Nx,1:Ny).*lap_Phi3(1:Nx,1:Ny)-Phi3(1:Nx,1:Ny).*lap_Phi2(1:Nx,1:Ny)))); 
                                

    termi33_2(1:Nx,1:Ny)=(16/(pi^2))*(sigma23*Phi2(1:Nx,1:Ny)+sigma13*Phi1(1:Nx,1:Ny)+sigma123*Phi1(1:Nx,1:Ny).*Phi2(1:Nx,1:Ny))/epsilon;

for i = 1 : Nx
    for j = 1 : Ny

        jp = j + 1;
        jm = j - 1;
        ip = i + 1;
        im = i - 1;

        if(im == 0)
           im=Nx;
        end
        if(ip == (Nx+1))
           ip=1;
        end
        if(jm == 0) 
           jm = 1;
        end
        if(jp == (Ny+1))
           jp=Ny;
        end
        % if(jm == 0) 
        %    jm = Ny;
        % end
        % if(jp == (Ny+1))
        %    jp=1;
        % end
        S1 = Phi1(ip, j) + Phi1(im, j) + Phi1(i, jp) + Phi1(i, jm);
        S2 = Phi2(ip, j) + Phi2(im, j) + Phi2(i, jp) + Phi2(i, jm);
        S3 = Phi3(ip, j) + Phi3(im, j) + Phi3(i, jp) + Phi3(i, jm);
        if(S1>tol && S1<4-tol)
              W1=1;
        else
              W1=0; 
        end
        if(S2>tol&& S2<4-tol)
            W2=1;
        else
            W2=0;   
        end
        if(S3>tol && S3<4-tol)
            W3=1;
        else
            W3=0;  
        end
         Lan_ta=W3+W2+W1;

        if ((Phi1(i, j) == 1 && Phi2(i, j) == 0 && Phi3(i, j) == 0) || (Phi1(i, j) == 0 && Phi2(i, j) == 1 && Phi3(i, j) == 0) || (Phi1(i, j) == 0 && Phi2(i, j) == 0 && Phi3(i, j) == 1) )

               
                if ((W1==1&&W2==1&&W3==0)||(W1==1&&W2==0&&W3==1))
                     Phi1_dot(i, j) = dt * (1.0 / epsilon) * (1.0 / Lan_ta)* ((1.0 / tau12)*W1*W2 * (termg(i, j) - termh(i, j)) + (1.0 / tau13) *W1*W3*(termg(i, j) - termi(i, j)));
     
                elseif(W1==1&&W2==1&&W3==1)
                 
                    Phi1_dot(i, j) = dt * (1.0 / epsilon)*(1.0 / Tau123)*(termg33(i,j)-termg33_2(i,j)-gr_po1(i,j)*dh_Phi1(i,j)-(N)*(termg33(i,j)-termg33_2(i,j)-gr_po1(i,j)*dh_Phi1(i,j)+termh33(i,j)-termh33_2(i,j)-gr_po2(i,j)*dh_Phi2(i,j)+termi33(i,j)-termi33_2(i,j)-gr_po3(i,j)*dh_Phi3(i,j)));
                else
                   Phi1_dot(i,j)=0; 
                end
                    

                


               


        elseif ( (Phi1(i, j) > 0 && Phi1(i, j) < 1) && (Phi2(i, j) > 0 && Phi2(i, j) < 1) && (Phi3(i, j) == 0) )
                  
            Phi1_dot(i, j) = dt * (1.0 / epsilon) * (1.0 / 2)* ((1.0 / tau12)*1*1 * (termg(i, j) - termh(i, j)) + (1.0 / tau13) *1*0*(termg(i, j) - termi(i, j)));%
                  
                  
        elseif ( (Phi1(i, j) > 0 && Phi1(i, j) < 1) && (Phi3(i, j) > 0 && Phi3(i, j) < 1) && (Phi2(i, j) == 0) )
                  
            Phi1_dot(i, j) = dt * (1.0 / epsilon) * (1.0 / 2)* ((1.0 / tau12)*1*0 * (termg(i, j) - termh(i, j)) + (1.0 / tau13) *1*1*(termg(i, j) - termi(i, j)));
                  

                  
        elseif ( (Phi1(i, j) > 0 && Phi1(i, j) < 1) && (Phi2(i, j) > 0 && Phi2(i, j) < 1) && (Phi3(i, j) > 0 && Phi3(i, j) < 1) )
           
            Phi1_dot(i, j) = dt * (1.0 / epsilon)*(1.0 / Tau123)*(termg33(i,j)-termg33_2(i,j)-gr_po1(i,j)*dh_Phi1(i,j)-(N)*(termg33(i,j)-termg33_2(i,j)-gr_po1(i,j)*dh_Phi1(i,j)+termh33(i,j)-termh33_2(i,j)-gr_po2(i,j)*dh_Phi2(i,j)+termi33(i,j)-termi33_2(i,j)-gr_po3(i,j)*dh_Phi3(i,j)));
        else
                   Phi1_dot(i, j) = 0;
       
        end



        if ((Phi1(i, j) == 1 && Phi2(i, j) == 0 && Phi3(i, j) == 0) || (Phi1(i, j) == 0 && Phi2(i, j) == 1 && Phi3(i, j) == 0) || (Phi1(i, j) == 0 && Phi2(i, j) == 0 && Phi3(i, j) == 1) )


            
                if ((W1==1&&W2==1&&W3==0)||(W1==0&&W2==1&&W3==1))

                    Phi2_dot(i, j) = dt * (1.0 / epsilon) * (1.0 /Lan_ta )* ((1.0 / tau12)*W1*W2 * (termh(i, j) - termg(i, j)) + (1.0 / tau23) *W3*W2*(termh(i, j) - termi(i, j)));
                  
                elseif(W1==1&&W2==1&&W3==1)
    
                Phi2_dot(i, j) = dt * (1.0 / epsilon)*(1.0 / Tau123)*(termh33(i,j)-termh33_2(i,j)-gr_po2(i,j)*dh_Phi2(i,j)-(N)*(termg33(i,j)-termg33_2(i,j)-gr_po1(i,j)*dh_Phi1(i,j)+termh33(i,j)-termh33_2(i,j)-gr_po2(i,j)*dh_Phi2(i,j)+termi33(i,j)-termi33_2(i,j)-gr_po3(i,j)*dh_Phi3(i,j)));
                else
                    Phi2_dot(i, j) =0;
                end
                

            

            
            



        elseif ( (Phi1(i, j) > 0 && Phi1(i, j) < 1) && (Phi2(i, j) > 0 && Phi2(i, j) < 1) && (Phi3(i, j) == 0) )

                  
                   Phi2_dot(i, j) = dt * (1.0 / epsilon) * (1.0 / 2)* ((1.0 / tau12)*1*1 * (termh(i, j) - termg(i, j)) + (1.0 / tau23) *0*1*(termh(i, j) - termi(i, j)));
                  

                  
        elseif ( (Phi2(i, j) > 0 && Phi2(i, j) < 1) && (Phi3(i, j) > 0 && Phi3(i, j) < 1) && (Phi1(i, j) == 0) )
                   
                   Phi2_dot(i, j) = dt * (1.0 / epsilon) * (1.0 / 2)* ((1.0 / tau12)*0*1 * (termh(i, j) - termg(i, j)) + (1.0 / tau23) *1*1*(termh(i, j) - termi(i, j)));



                    
        elseif ( (Phi1(i, j) > 0 && Phi1(i, j) < 1) && (Phi2(i, j) > 0 && Phi2(i, j) < 1) && (Phi3(i, j) > 0 && Phi3(i, j) < 1) )
                    % W1=1; W2=1;W3=1;
                    % tau123=(tau12*Phi1(i, j)*Phi2(i, j)+tau13*Phi1(i, j)*Phi3(i, j)+tau23*Phi2(i, j)*Phi3(i, j))/(Phi1(i, j)*Phi2(i, j)+Phi1(i, j)*Phi3(i, j)+Phi2(i, j)*Phi3(i, j));
                    % Phi2_dot(i, j) = dt * (1.0 / epsilon) * (1.0 / 3)* ((1.0 / tau12)*W1*W2 * (termh(i, j) - termg(i, j)) + (1.0 / tau23) *W3*W2*(termh(i, j) - termi(i, j)));
                    Phi2_dot(i, j) = dt * (1.0 / epsilon)*(1.0 / Tau123)*(termh33(i,j)-termh33_2(i,j)-gr_po2(i,j)*dh_Phi2(i,j)-(N)*(termg33(i,j)-termg33_2(i,j)-gr_po1(i,j)*dh_Phi1(i,j)+termh33(i,j)-termh33_2(i,j)-gr_po2(i,j)*dh_Phi2(i,j)+termi33(i,j)-termi33_2(i,j)-gr_po3(i,j)*dh_Phi3(i,j)));
         else
                Phi2_dot(i, j) = 0;
       
        end



        if ((Phi1(i, j) == 1 && Phi2(i, j) == 0 && Phi3(i, j) == 0) || (Phi1(i, j) == 0 && Phi2(i, j) == 1 && Phi3(i, j) == 0) || (Phi1(i, j) == 0 && Phi2(i, j) == 0 && Phi3(i, j) == 1) )


            
                    if ((W1==1&&W2==0&&W3==1)||(W1==0&&W2==1&&W3==1))
                    
                    Phi3_dot(i, j) = dt * (1.0 / epsilon) * (1.0 / Lan_ta)* ((1.0 / tau13)*W1*W3 * (termi(i, j) - termg(i, j)) + (1.0 / tau23) *W3*W2*(termi(i, j) - termh(i, j)));
                    
                    elseif(W1==1&&W2==1&&W3==1)
                    
                    % Phi3_dot(i, j) = dt * (1.0 / epsilon)*(1.0 / Tau123)*(termi33(i,j)-gr_po3(i,j)*dh_Phi3(i,j)-(1.0/3)*(termg33(i,j)-gr_po1(i,j)*dh_Phi1(i,j)+termh33(i,j)-gr_po2(i,j)*dh_Phi2(i,j)+termi33(i,j)-gr_po3(i,j)*dh_Phi3(i,j)));
                    Phi3_dot(i, j) =-Phi1_dot(i, j) -Phi2_dot(i, j) ;
                    else
                       Phi3_dot(i, j) =0; 
                    
                    
                    end
            
            



    elseif ( (Phi1(i, j) > 0 && Phi1(i, j) < 1) && (Phi3(i, j) > 0 && Phi3(i, j) < 1) && (Phi2(i, j) == 0) )

             
                    Phi3_dot(i, j) = dt * (1.0 / epsilon) * (1.0 / 2)* ((1.0 / tau13)*1*1 * (termi(i, j) - termg(i, j)) + (1.0 / tau23) *1*0*(termi(i, j) - termh(i, j)));
            

                   
                    
           

      elseif ( (Phi2(i, j) > 0 && Phi2(i, j) < 1) && (Phi3(i, j) > 0 && Phi3(i, j) < 1) && (Phi1(i, j) == 0) )

               
                    Phi3_dot(i, j) = dt * (1.0 / epsilon) * (1.0 / 2)* ((1.0 / tau13)*0*1 * (termi(i, j) - termg(i, j)) + (1.0 / tau23) *1*1*(termi(i, j) - termh(i, j)));
                

                    % Phi3_dot(i, j) = dt * (1.0 / epsilon)*(1.0 / Tau123)*(termi33(i,j)-gr_po3(i,j)*dh_Phi3(i,j)-(1.0/3)*(termg33(i,j)-gr_po1(i,j)*dh_Phi1(i,j)+termh33(i,j)-gr_po2(i,j)*dh_Phi2(i,j)+termi33(i,j)-gr_po3(i,j)*dh_Phi3(i,j)));
                 

        elseif ( (Phi1(i, j) > 0 && Phi1(i, j) < 1) && (Phi2(i, j) > 0 && Phi2(i, j) < 1) && (Phi3(i, j) > 0 && Phi3(i, j) < 1) )

            
                    tau123=(tau12*Phi1(i, j)*Phi2(i, j)+tau13*Phi1(i, j)*Phi3(i, j)+tau23*Phi2(i, j)*Phi3(i, j))/(Phi1(i, j)*Phi2(i, j)+Phi1(i, j)*Phi3(i, j)+Phi2(i, j)*Phi3(i, j));
                    % Phi3_dot(i, j) = dt * (1.0 / epsilon) * (1.0 / 3)* ((1.0 / tau13)*W1*W3 * (termi3(i, j) - termg1(i, j)) + (1.0 / tau23) *W3*W2*(termi3(i, j) - termh2(i, j)));
                    % Phi3_dot(i, j) = dt * (1.0 / epsilon)*(1.0 / tau123)*(termi33(i,j)-gr_po3(i,j)*dh_Phi3(i,j)-(1.0/3)*(termg33(i,j)-gr_po1(i,j)*dh_Phi1(i,j)+termh33(i,j)-gr_po2(i,j)*dh_Phi2(i,j)+termi33(i,j)-gr_po3(i,j)*dh_Phi3(i,j)));
                    Phi3_dot(i, j) =-Phi1_dot(i, j) -Phi2_dot(i, j) ;
         else

            Phi3_dot(i, j) = 0;

        
        end
        Phi_dot(i,j)=Phi1_dot(i, j)+Phi2_dot(i, j)+Phi3_dot(i, j);


     
        
    end
    end


for i = 1 : Nx
    for j = 1 : Ny
        
        Phi1(i, j) = Phi1(i, j) + Phi1_dot(i, j);
        Phi2(i, j) = Phi2(i, j) + Phi2_dot(i, j);
        
        Phi3(i, j) = Phi3(i, j) + Phi3_dot(i, j);
       
        
        if (Phi1(i, j) >=1-tol)
            Phi1(i, j) = 1.0;  end

        if (Phi2(i, j) >= 1-tol)
            Phi2(i, j) = 1.0;  end

        if (Phi3(i, j) >= 1-tol)
            Phi3(i, j) = 1.0;  end
        
        if (Phi1(i, j) <= tol)
            Phi1(i, j) = 0.0;  end
        
        if (Phi2(i, j) <= tol)
            Phi2(i, j) = 0.0;  end
        
        if (Phi3(i, j) <= tol)
            Phi3(i, j) = 0.0;  end
        
          Sum_phi = Phi3(i, j)+  Phi2(i, j) +  Phi1(i, j);
          Phi1(i, j)=Phi1(i, j)/Sum_phi;
          Phi2(i, j)=Phi2(i, j)/Sum_phi;
          Phi3(i, j)=Phi3(i, j)/Sum_phi;
          Phi(i,j) = Phi1(i, j) + Phi2(i, j)+ Phi3(i, j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end

for i = 1 : Nx
    for j = 1 : Ny

        

        
     
    if((Phi1(i, j) > 0 && Phi1(i, j) < 1) && (Phi2(i, j) > 0 && Phi2(i, j) < 1) && (Phi3(i, j) > 0 && Phi3(i, j) < 1))
         dh_Phi1(i,j)=6*Phi1(i,j)*(1-Phi1(i,j))+2*Phi2(i,j)*Phi3(i,j);
         dh_Phi2(i,j)=6*Phi2(i,j)*(1-Phi2(i,j))+2*Phi1(i,j)*Phi3(i,j); 
         dh_Phi3(i,j)=6*Phi3(i,j)*(1-Phi3(i,j))+2*Phi1(i,j)*Phi2(i,j);   
    else
    dh_Phi1(i,j)=6*Phi1(i,j)*(1-Phi1(i,j));
    dh_Phi2(i,j)=6*Phi2(i,j)*(1-Phi2(i,j)); 
    dh_Phi3(i,j)=6*Phi3(i,j)*(1-Phi3(i,j));
        
     end
    end
    end



        
        hp_Phi1(1:Nx,1:Ny)=(Phi1(1:Nx,1:Ny).^2).*(3-2*Phi1(1:Nx,1:Ny))+2*Phi1(1:Nx,1:Ny).*Phi2(1:Nx,1:Ny).*Phi3(1:Nx,1:Ny);
        hp_Phi2(1:Nx,1:Ny)=(Phi2(1:Nx,1:Ny).^2).*(3-2*Phi2(1:Nx,1:Ny))+2*Phi1(1:Nx,1:Ny).*Phi2(1:Nx,1:Ny).*Phi3(1:Nx,1:Ny);
        hp_Phi3(1:Nx,1:Ny)=(Phi3(1:Nx,1:Ny).^2).*(3-2*Phi3(1:Nx,1:Ny))+2*Phi1(1:Nx,1:Ny).*Phi2(1:Nx,1:Ny).*Phi3(1:Nx,1:Ny);

        % hp_Phi1(1:Nx,1:Ny)=(Phi1(1:Nx,1:Ny).^2).*(3-2*Phi1(1:Nx,1:Ny));
        % hp_Phi2(1:Nx,1:Ny)=(Phi2(1:Nx,1:Ny).^2).*(3-2*Phi2(1:Nx,1:Ny));
        % hp_Phi3(1:Nx,1:Ny)=(Phi3(1:Nx,1:Ny).^2).*(3-2*Phi3(1:Nx,1:Ny));

        % HH1(1:Nx,1:Ny)=hhp_Phi1(1:Nx,1:Ny)-hp_Phi1(1:Nx,1:Ny);
        % HH2(1:Nx,1:Ny)=hhp_Phi2(1:Nx,1:Ny)-hp_Phi2(1:Nx,1:Ny);
        % HH3(1:Nx,1:Ny)=hhp_Phi3(1:Nx,1:Ny)-hp_Phi3(1:Nx,1:Ny);
        % 

       

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Update of Diffusion Potential 


        M11_Phi(1:Nx,1:Ny)= D11*dCom1_dmu1*hp_Phi1(1:Nx,1:Ny)+D21*dCom2_dmu1*hp_Phi2(1:Nx,1:Ny)...
                          + D31*dCom3_dmu1*hp_Phi3(1:Nx,1:Ny);  
        M12_Phi(1:Nx,1:Ny)= D12*dCom1_dmu2*hp_Phi1(1:Nx,1:Ny)+D22*dCom2_dmu2*hp_Phi2(1:Nx,1:Ny)...
                          + D32*dCom3_dmu2*hp_Phi3(1:Nx,1:Ny); 
        M21_Phi(1:Nx,1:Ny)= D13*dCom1_dmu3*hp_Phi1(1:Nx,1:Ny)+D23*dCom2_dmu3*hp_Phi2(1:Nx,1:Ny)...
                          + D33*dCom3_dmu3*hp_Phi3(1:Nx,1:Ny); 
        M22_Phi(1:Nx,1:Ny)= D14*dCom1_dmu4*hp_Phi1(1:Nx,1:Ny)+D24*dCom2_dmu4*hp_Phi2(1:Nx,1:Ny)...
                          + D34*dCom3_dmu4*hp_Phi3(1:Nx,1:Ny); 

for i=1:Nx 
    for j=1:Ny
         jp=j+1;
         jm=j-1;
         ip=i+1;
         im=i-1;
        if(im == 0)
           im=Nx;
        end
        if(ip == (Nx+1))
           ip=1;
        end
        if(jm == 0) 
           jm =1;
        end
        if(jp == (Ny+1))
           jp=Ny;
        end
        % if(jm == 0) 
        %    jm = Ny;
        % end
        % if(jp == (Ny+1))
        %    jp=1;
        % end
        

  fenzi1(i,j)=(((M11_Phi(ip,j)+M11_Phi(i,j))*(miu(ip,j)-miu(i,j))...
                -(M11_Phi(i,j)+M11_Phi(im,j))*(miu(i,j)-miu(im,j)))/(2*(dx^2)))+...
              (((M11_Phi(i,jp)+M11_Phi(i,j))*(miu(i,jp)-miu(i,j))...
                -(M11_Phi(i,j)+M11_Phi(i,jm))*(miu(i,j)-miu(i,jm)))/(2*(dy^2)))+...
              (((M12_Phi(ip,j)+M12_Phi(i,j))*(miu2(ip,j)-miu2(i,j))...
                -(M12_Phi(i,j)+M12_Phi(im,j))*(miu2(i,j)-miu2(im,j)))/(2*(dx^2)))+...
              (((M12_Phi(i,jp)+M12_Phi(i,j))*(miu2(i,jp)-miu2(i,j))...
                -(M12_Phi(i,j)+M12_Phi(i,jm))*(miu2(i,j)-miu2(i,jm)))/(2*(dy^2)));
              
  fenzi12(i,j)=(((M21_Phi(ip,j)+M21_Phi(i,j))*(miu(ip,j)-miu(i,j))...
                -(M21_Phi(i,j)+M21_Phi(im,j))*(miu(i,j)-miu(im,j)))/(2*(dx^2)))+...
               (((M21_Phi(i,jp)+M21_Phi(i,j))*(miu(i,jp)-miu(i,j))...
                -(M21_Phi(i,j)+M21_Phi(i,jm))*(miu(i,j)-miu(i,jm)))/(2*(dy^2)))+...
               (((M22_Phi(ip,j)+M22_Phi(i,j))*(miu2(ip,j)-miu2(i,j))...
                -(M22_Phi(i,j)+M22_Phi(im,j))*(miu2(i,j)-miu2(im,j)))/(2*(dx^2)))+...
               (((M22_Phi(i,jp)+M22_Phi(i,j))*(miu2(i,jp)-miu2(i,j))...
                -(M22_Phi(i,j)+M22_Phi(i,jm))*(miu2(i,j)-miu2(i,jm)))/(2*(dy^2)));
              
    end
end  
%%%涉及到插值函数求导，dh_Phi1 等求导之后的函数，使得三相点外还起作用
%%%（本来h_Phi1中2*Phi1(1:Nx,1:Ny).*Phi2(1:Nx,1:Ny).*Phi3(1:Nx,1:Ny))限制其只在三相点处起作用，但是一求导，限制消失了



        fenzi2(1 : Nx, 1 : Ny) = Com11(1 : Nx, 1 : Ny) .* dh_Phi1(1 : Nx, 1 : Ny) .* Phi1_dot(1 : Nx, 1 : Ny)/dt ...
                               + Com21(1 : Nx, 1 : Ny) .* dh_Phi2(1 : Nx, 1 : Ny) .* Phi2_dot(1 : Nx, 1 : Ny)/dt ...
                               + Com31(1 : Nx, 1 : Ny) .* dh_Phi3(1 : Nx, 1 : Ny) .* Phi3_dot(1 : Nx, 1 : Ny)/dt ;
         
        fenzi22(1 : Nx, 1 : Ny) = Com12(1 : Nx, 1 : Ny) .* dh_Phi1(1 : Nx, 1 : Ny) .* Phi1_dot(1 : Nx, 1 : Ny)/dt ...
                                + Com22(1 : Nx, 1 : Ny) .* dh_Phi2(1 : Nx, 1 : Ny) .* Phi2_dot(1 : Nx, 1 : Ny)/dt ...
                                + Com32(1 : Nx, 1 : Ny) .* dh_Phi3(1 : Nx, 1 : Ny) .* Phi3_dot(1 : Nx, 1 : Ny)/dt ;
     
        fenmu(1 : Nx, 1 : Ny) = (dCom1_dmu1+dCom1_dmu2)*hp_Phi1(1:Nx,1:Ny)+...
                                (dCom2_dmu1+dCom2_dmu2)*hp_Phi2(1:Nx,1:Ny)+...
                                (dCom3_dmu1+dCom3_dmu2)*hp_Phi3(1:Nx,1:Ny);
        fenmu1(1 : Nx, 1 : Ny) = (dCom1_dmu3+dCom1_dmu4)*hp_Phi1(1:Nx,1:Ny)+...
                                (dCom2_dmu3+dCom2_dmu4)*hp_Phi2(1:Nx,1:Ny)+...
                                (dCom3_dmu4+dCom3_dmu3)*hp_Phi3(1:Nx,1:Ny);

     
     a(1:Nx,1:Ny)=(hp_Phi1(1:Nx,1:Ny)*dCom1_dmu1+hp_Phi2(1:Nx,1:Ny)*dCom2_dmu1+hp_Phi3(1:Nx,1:Ny)*dCom3_dmu1);
     b(1:Nx,1:Ny)=(hp_Phi1(1:Nx,1:Ny)*dCom1_dmu2+hp_Phi2(1:Nx,1:Ny)*dCom2_dmu2+hp_Phi3(1:Nx,1:Ny)*dCom3_dmu2);
     c(1:Nx,1:Ny)=(hp_Phi1(1:Nx,1:Ny)*dCom1_dmu3+hp_Phi2(1:Nx,1:Ny)*dCom2_dmu3+hp_Phi3(1:Nx,1:Ny)*dCom3_dmu3);
     d(1:Nx,1:Ny)=(hp_Phi1(1:Nx,1:Ny)*dCom1_dmu4+hp_Phi2(1:Nx,1:Ny)*dCom2_dmu4+hp_Phi3(1:Nx,1:Ny)*dCom3_dmu4);


    
%%%2025.12.14
for i=1:Nx 
    for j=1:Ny
         
        if(Phi1(i,j)==1||Phi1(i,j)==0||Phi3(i,j)==1||Phi3(i,j)==0)
            J11_1(i,j)=0;
                 J21_1(i,j)=0;
        elseif((Phi1(i,j)>0&&Phi1(i,j)<1)&&(Phi3(i,j)>0&&Phi3(i,j)<1)&&((grad_x_Phi1(i,j)^2+grad_y_Phi1(i,j)^2)>0)&&((grad_x_Phi3(i,j)^2+grad_y_Phi3(i,j)^2)>0))
            J11_1(i,j)=((pi*epsilon)/4)*(((hp_Phi1(i,j)*(1-hp_Phi1(i,j)))/sqrt(Phi1(i,j)*(1-Phi1(i,j)))))* Phi1_dot(i,j)./dt.*(((g_Phi13(i,j))/((sqrt(g_Phi11(i,j))*sqrt(g_Phi33(i,j))))))*(Com31(i,j)-Com11(i,j));
            J21_1(i,j)=((pi*epsilon)/4)*(((hp_Phi1(i,j)*(1-hp_Phi1(i,j)))/sqrt(Phi1(i,j)*(1-Phi1(i,j)))))* Phi1_dot(i,j)./dt.*(((g_Phi13(i,j))/((sqrt(g_Phi11(i,j))*sqrt(g_Phi33(i,j))))))*(Com32(i,j)-Com12(i,j));
        end
        if(Phi2(i,j)==1||Phi2(i,j)==0||Phi3(i,j)==1||Phi3(i,j)==0)
         J12_1(i,j)=0;
                 J22_1(i,j)=0; 
        elseif((Phi2(i,j)>0&&Phi2(i,j)<1)&&(Phi3(i,j)>0&&Phi3(i,j)<1)&&((grad_x_Phi2(i,j)^2+grad_y_Phi2(i,j)^2)>0)&&((grad_x_Phi3(i,j)^2+grad_y_Phi3(i,j)^2)>0))
                 J12_1(i,j)=((pi*epsilon)/4)*(((hp_Phi2(i,j)*(1-hp_Phi2(i,j)))/sqrt(Phi2(i,j)*(1-Phi2(i,j)))))* Phi2_dot(i,j)./dt.*(((g_Phi23(i,j))/((sqrt(g_Phi22(i,j))*sqrt(g_Phi33(i,j))))))*(Com31(i,j)-Com21(i,j));
                 J22_1(i,j)=((pi*epsilon)/4)*(((hp_Phi2(i,j)*(1-hp_Phi2(i,j)))/sqrt(Phi2(i,j)*(1-Phi2(i,j)))))* Phi2_dot(i,j)./dt.*(((g_Phi23(i,j))/((sqrt(g_Phi22(i,j)))*sqrt(g_Phi33(i,j)))))*(Com32(i,j)-Com22(i,j));
        end
    end
 end
%%%%end
     
for i=1:Nx 
    for j=1:Ny
         jp=j+1;
         jm=j-1;
         ip=i+1;
         im=i-1;
        if(im == 0)
           im=Nx;
        end
        if(ip == (Nx+1))
           ip=1;
        end
        if(jm == 0) 
           jm = 1;
        end
        if(jp == (Ny+1))
           jp=Ny;
        end
        % if(jm == 0) 
        %    jm = Ny;
        % end
        % if(jp == (Ny+1))
        %    jp=1;
        % end
             if(Phi1(i,j)==1||Phi1(i,j)==0||Phi3(i,j)==1||Phi3(i,j)==0)
                 J11_1(i,j)=0;
                 J21_1(i,j)=0;
                 % J11(i,j)=0 ; 
                 % J21(i,j)=0 ;
                 x2_J11(i,j)=0;
                 x1_J11(i,j)=0;
                 y2_J11(i,j)=0;
                 y1_J11(i,j)=0;
                 x2_J21(i,j)=0;
                 x1_J21(i,j)=0;
                 y2_J21(i,j)=0;
                 y1_J21(i,j)=0;
             elseif((Phi1(i,j)>0&&Phi1(i,j)<1)&&(Phi3(i,j)>0&&Phi3(i,j)<1)&&((grad_x_Phi1(i,j)^2+grad_y_Phi1(i,j)^2)>0)&&((grad_x_Phi3(i,j)^2+grad_y_Phi3(i,j)^2)>0))
                %J11(i+1/2,j)
                
                 x2_J11_1(i,j)=(J11_1(i,j)+J11_1(ip,j))/2;
                numer_x2_J11(i,j)=(Phi1(ip,j)-Phi1(i,j))/dx;

                denom_x2_J11e_x(i,j)=(Phi1(ip,j)-Phi1(i,j))/dx;
                denom_x2_J11e_y(i,j)=((Phi1(i,jp)-Phi1(i,jm))/(2*dy)+(Phi1(ip,jp)-Phi1(ip,jm))/(2*dy))/2;

                denom_x2_J11(i,j)=sqrt(denom_x2_J11e_x(i,j)^2+denom_x2_J11e_y(i,j)^2);

                if(denom_x2_J11(i,j)==0)
                    denom_x2_J11(i,j)=tol;
                else
                    denom_x2_J11(i,j)=sqrt(denom_x2_J11e_x(i,j)^2+denom_x2_J11e_y(i,j)^2);
                end
                
                x2_J11(i,j)=x2_J11_1(i,j)*numer_x2_J11(i,j)/denom_x2_J11(i,j);
                %J11(i-1/2,j)

                x1_J11_1(i,j)=(J11_1(i,j)+J11_1(im,j))/2;
                numer_x1_J11(i,j)=(Phi1(i,j)-Phi1(im,j))/dx;

                denom_x1_J11w_x(i,j)=(Phi1(i,j)-Phi1(im,j))/dx;
                denom_x1_J11w_y(i,j)=((Phi1(i,jp)-Phi1(i,jm))/(2*dy)+(Phi1(im,jp)-Phi1(im,jm))/(2*dy))/2;

                denom_x1_J11(i,j)=sqrt(denom_x1_J11w_x(i,j)^2+denom_x1_J11w_y(i,j)^2);
                if(denom_x1_J11(i,j)==0)
                    denom_x1_J11(i,j)=tol;
                else
                    denom_x1_J11(i,j)=sqrt(denom_x1_J11w_x(i,j)^2+denom_x1_J11w_y(i,j)^2);
                end
                
                x1_J11(i,j)=x1_J11_1(i,j)*numer_x1_J11(i,j)/denom_x1_J11(i,j);

                %J11(i,j+1/2)
                y2_J11_1(i,j)=(J11_1(i,j)+J11_1(i,jp))/2;
                numer_y2_J11(i,j)=(Phi1(i,jp)-Phi1(i,j))/dy;


                denom_y2_J11n_x(i,j)=((Phi1(ip,j)-Phi1(im,j))/(2*dx)+(Phi1(ip,jp)-Phi1(im,jp))/(2*dx))/2;
                denom_y2_J11n_y(i,j)=(Phi1(i,jp)-Phi1(i,j))/dy;

                denom_y2_J11(i,j)=sqrt(denom_y2_J11n_x(i,j)^2+denom_y2_J11n_y(i,j)^2);

                if(denom_y2_J11(i,j)==0)
                   denom_y2_J11(i,j)=tol;
                else
                    denom_y2_J11(i,j)=sqrt(denom_y2_J11n_x(i,j)^2+denom_y2_J11n_y(i,j)^2);
                end
                
                y2_J11(i,j)=y2_J11_1(i,j)*numer_y2_J11(i,j)/denom_y2_J11(i,j);

                %J11(i,j-1/2);

                y1_J11_1(i,j)=(J11_1(i,j)+J11_1(i,jm))/2;
                numer_y1_J11(i,j)=(Phi1(i,j)-Phi1(i,jm))/dy;

                denom_y1_J11s_x(i,j)=((Phi1(ip,j)-Phi1(im,j))/(2*dx)+(Phi1(ip,jm)-Phi1(im,jm))/(2*dx))/2;
                denom_y1_J11s_y(i,j)=(Phi1(i,j)-Phi1(i,jm))/dy;

                denom_y1_J11(i,j)=sqrt(denom_y1_J11s_x(i,j)^2+denom_y1_J11s_y(i,j)^2);

                if(denom_y1_J11(i,j)==0)
                    denom_y1_J11(i,j)=tol;
                else
                    denom_y1_J11(i,j)=sqrt(denom_y1_J11s_x(i,j)^2+denom_y1_J11s_y(i,j)^2);
                end
                
                y1_J11(i,j)=y1_J11_1(i,j)*numer_y1_J11(i,j)/denom_y1_J11(i,j);




                

                
               
               
                
             
                
                %J21(i+1/2,j)
                x2_J21_1(i,j)=(J21_1(i,j)+J21_1(ip,j))/2;
                numer_x2_J21(i,j)=(Phi1(ip,j)-Phi1(i,j))/dx;

                denom_x2_J21e_x(i,j)=(Phi1(ip,j)-Phi1(i,j))/dx;
                denom_x2_J21e_y(i,j)=((Phi1(i,jp)-Phi1(i,jm))/(2*dy)+(Phi1(ip,jp)-Phi1(ip,jm))/(2*dy))/2;

                denom_x2_J21(i,j)=sqrt(denom_x2_J21e_x(i,j)^2+denom_x2_J21e_y(i,j)^2);

                if(denom_x2_J21(i,j)==0)
                    denom_x2_J21(i,j)=tol;
                else
                    denom_x2_J21(i,j)=sqrt(denom_x2_J21e_x(i,j)^2+denom_x2_J21e_y(i,j)^2);
                end
                
                x2_J21(i,j)=x2_J21_1(i,j)*numer_x2_J21(i,j)/denom_x2_J21(i,j);
                
                %J21(i-1/2,j)
                x1_J21_1(i,j)=(J21_1(i,j)+J21_1(im,j))/2;
                numer_x1_J21(i,j)=(Phi1(i,j)-Phi1(im,j))/dx;

                denom_x1_J21w_x(i,j)=(Phi1(i,j)-Phi1(im,j))/dx;
                denom_x1_J21w_y(i,j)=((Phi1(i,jp)-Phi1(i,jm))/(2*dy)+(Phi1(im,jp)-Phi1(im,jm))/(2*dy))/2;

                denom_x1_J21(i,j)=sqrt(denom_x1_J21w_x(i,j)^2+denom_x1_J21w_y(i,j)^2);
                
                if(denom_x1_J21(i,j)==0)
                    denom_x1_J21(i,j)=tol;
                else
                    denom_x1_J21(i,j)=sqrt(denom_x1_J21w_x(i,j)^2+denom_x1_J21w_y(i,j)^2);
                end
                
                x1_J21(i,j)=x1_J21_1(i,j)*numer_x1_J21(i,j)/denom_x1_J21(i,j);
                
                %J21(i,j+1/2)
                y2_J21_1(i,j)=(J21_1(i,j)+J21_1(i,jp))/2;
                numer_y2_J21(i,j)=(Phi1(i,jp)-Phi1(i,j))/dy;

                denom_y2_J21n_x(i,j)=((Phi1(ip,j)-Phi1(im,j))/(2*dx)+(Phi1(ip,jp)-Phi1(im,jp))/(2*dx))/2;
                denom_y2_J21n_y(i,j)=(Phi1(i,jp)-Phi1(i,j))/dy;

                denom_y2_J21(i,j)=sqrt(denom_y2_J21n_x(i,j)^2+denom_y2_J21n_y(i,j)^2);

                 if(denom_y2_J21(i,j)==0)
                    denom_y2_J21(i,j)=tol;
                else
                    denom_y2_J21(i,j)=sqrt(denom_y2_J21n_x(i,j)^2+denom_y2_J21n_y(i,j)^2);
                end
                
                y2_J21(i,j)=y2_J21_1(i,j)*numer_y2_J21(i,j)/denom_y2_J21(i,j);

                %J21(i,j-1/2);

                y1_J21_1(i,j)=(J21_1(i,j)+J21_1(i,jm))/2;
                numer_y1_J21(i,j)=(Phi1(i,j)-Phi1(i,jm))/dy;

                denom_y1_J21s_x(i,j)=((Phi1(ip,j)-Phi1(im,j))/(2*dx)+(Phi1(ip,jm)-Phi1(im,jm))/(2*dx))/2;
                denom_y1_J21s_y(i,j)=(Phi1(i,j)-Phi1(i,jm))/dy;

                denom_y1_J21(i,j)=sqrt(denom_y1_J21s_x(i,j)^2+denom_y1_J21s_y(i,j)^2);

                if(denom_y1_J21(i,j)==0)
                    denom_y1_J21(i,j)=tol;
                else
                    denom_y1_J21(i,j)=sqrt(denom_y1_J21s_x(i,j)^2+denom_y1_J21s_y(i,j)^2);
                end
                
                y1_J21(i,j)=y1_J21_1(i,j)*numer_y1_J21(i,j)/denom_y1_J21(i,j);





             end
             if(Phi2(i,j)==1||Phi2(i,j)==0||Phi3(i,j)==1||Phi3(i,j)==0)
                 J12_1(i,j)=0;
                 J22_1(i,j)=0;
                 % J12(i,j)=0;
                 % J22(i,j)=0;
                 x2_J12(i,j)=0;
                 x1_J12(i,j)=0;
                 y2_J12(i,j)=0;
                 y1_J12(i,j)=0;
                 x2_J22(i,j)=0;
                 x1_J22(i,j)=0;
                 y2_J22(i,j)=0;
                 y1_J22(i,j)=0;
                 

             elseif((Phi2(i,j)>0&&Phi2(i,j)<1)&&(Phi3(i,j)>0&&Phi3(i,j)<1)&&((grad_x_Phi2(i,j)^2+grad_y_Phi2(i,j)^2)>0)&&((grad_x_Phi3(i,j)^2+grad_y_Phi3(i,j)^2)>0))
                
                 %J12(i+1/2,j)

                x2_J12_1(i,j)=(J12_1(i,j)+J12_1(ip,j))/2;
                numer_x2_J12(i,j)=(Phi2(ip,j)-Phi2(i,j))/dx;

                denom_x2_J12e_x(i,j)=(Phi2(ip,j)-Phi2(i,j))/dx;
                denom_x2_J12e_y(i,j)=((Phi2(i,jp)-Phi2(i,jm))/(2*dy)+(Phi2(ip,jp)-Phi2(ip,jm))/(2*dy))/2;

                denom_x2_J12(i,j)=sqrt(denom_x2_J12e_x(i,j)^2+denom_x2_J12e_y(i,j)^2);

                if(denom_x2_J12(i,j)==0)
                    denom_x2_J12(i,j)=tol;
                else
                    denom_x2_J12(i,j)=sqrt(denom_x2_J12e_x(i,j)^2+denom_x2_J12e_y(i,j)^2);
                end
                
                x2_J12(i,j)=x2_J12_1(i,j)*numer_x2_J12(i,j)/denom_x2_J12(i,j);

                %J12(i-1/2,j)

                x1_J12_1(i,j)=(J12_1(i,j)+J12_1(im,j))/2;
                numer_x1_J12(i,j)=(Phi2(i,j)-Phi2(im,j))/dx;

                denom_x1_J12w_x(i,j)=((Phi2(i,j)-Phi2(im,j))/dx);
                denom_x1_J12w_y(i,j)=(((Phi2(i,jp)-Phi2(i,jm))/(2*dy)+(Phi2(im,jp)-Phi2(im,jm))/(2*dy))/2);

                denom_x1_J12(i,j)=sqrt(denom_x1_J12w_x(i,j)^2+denom_x1_J12w_y(i,j)^2);

                if(denom_x1_J12(i,j)==0)
                    denom_x1_J12(i,j)=tol;
                else
                    denom_x1_J12(i,j)=sqrt(denom_x1_J12w_x(i,j)^2+denom_x1_J12w_y(i,j)^2);
                end

                
                x1_J12(i,j)=x1_J12_1(i,j)*numer_x1_J12(i,j)/denom_x1_J12(i,j);

                %J12(i,j+1/2)
                
                y2_J12_1(i,j)=(J12_1(i,j)+J12_1(i,jp))/2;
                numer_y2_J12(i,j)=(Phi2(i,jp)-Phi2(i,j))/dy;

                denom_y2_J12n_x(i,j)=((Phi2(ip,j)-Phi2(im,j))/(2*dx)+(Phi2(ip,jp)-Phi2(im,jp))/(2*dx))/2;
                denom_y2_J12n_y(i,j)=(Phi2(i,jp)-Phi2(i,j))/dy;

                denom_y2_J12(i,j)=sqrt(denom_y2_J12n_x(i,j)^2+enom_y2_J12n_y(i,j)^2);

                if(denom_y2_J12(i,j)==0)
                    denom_y2_J12(i,j)=tol;
                else
                    denom_y2_J12(i,j)=sqrt(denom_y2_J12n_x(i,j)^2+enom_y2_J12n_y(i,j)^2);
                end

                
                y2_J12(i,j)=y2_J12_1(i,j)*numer_y2_J12(i,j)/denom_y2_J12(i,j);

                %J12(i,j-1/2);

                y1_J12_1(i,j)=(J12_1(i,j)+J12_1(i,jm))/2;
                numer_y1_J12(i,j)=(Phi2(i,j)-Phi2(i,jm))/dy;

                denom_y1_J12s_x(i,j)=((Phi2(ip,j)-Phi2(im,j))/(2*dx)+(Phi2(ip,jm)-Phi2(im,jm))/(2*dx))/2;
                denom_y1_J12s_y(i,j)=((Phi2(i,j)-Phi2(i,jm))/dy);

                denom_y1_J12(i,j)=sqrt(denom_y1_J12s_x(i,j)^2+denom_y1_J12s_y(i,j)^2);

                if(denom_y1_J12(i,j)==0)
                    denom_y1_J12(i,j)=tol;
                else
                    denom_y1_J12(i,j)=sqrt(denom_y1_J12s_x(i,j)^2+denom_y1_J12s_y(i,j)^2);
                end
                
                y1_J12(i,j)=y1_J12_1(i,j)*numer_y1_J12(i,j)/denom_y1_J12(i,j);




                % J12(i,j)=((pi*epsilon)/4)*(((hp_Phi2(i,j)*(1-hp_Phi2(i,j)))/sqrt(Phi2(i,j)*(1-Phi2(i,j)))))* Phi2_dot(i,j)./dt*(((g_Phi23(i,j))/((sqrt(g_Phi22(i,j))*sqrt(g_Phi33(i,j))))))*...
                %             (((Com31(i,j)-Com21(i,j))*grad_x_Phi2(i,j)+(Com31(i,j)-Com21(i,j))*grad_y_Phi2(i,j)))*(1/sqrt(g_Phi22(i,j)));

                J22_1(i,j)=((pi*epsilon)/4)*(((hp_Phi2(i,j)*(1-hp_Phi2(i,j)))/sqrt(Phi2(i,j)*(1-Phi2(i,j)))))* Phi2_dot(i,j)./dt*(((g_Phi23(i,j))/((sqrt(g_Phi22(i,j)))*sqrt(g_Phi33(i,j)))))*(Com32(i,j)-Com22(i,j));

                %J22(i+1/2,j)

                x2_J22_1(i,j)=(J22_1(i,j)+J22_1(ip,j))/2;
                numer_x2_J22(i,j)=(Phi2(ip,j)-Phi2(i,j))/dx;

                denom_x2_J22e_x(i,j)=((Phi2(ip,j)-Phi2(i,j))/dx);
                denom_x2_J22e_y(i,j)=((Phi2(i,jp)-Phi2(i,jm))/(2*dy)+(Phi2(ip,jp)-Phi2(ip,jm))/(2*dy))/2;

                denom_x2_J22(i,j)=sqrt(denom_x2_J22e_x(i,j)^2+denom_x2_J22e_y(i,j)^2);

                if(denom_x2_J22(i,j)==0)
                    denom_x2_J22(i,j)=tol;
                else
                    denom_x2_J22(i,j)=sqrt(denom_x2_J22e_x(i,j)^2+denom_x2_J22e_y(i,j)^2);
                end
                
                x2_J22(i,j)=x2_J22_1(i,j)*numer_x2_J22(i,j)/denom_x2_J22(i,j);

                %J22(i-1/2,j)

                x1_J22_1(i,j)=(J22_1(i,j)+J22_1(im,j))/2;
                numer_x1_J22(i,j)=(Phi2(i,j)-Phi2(im,j))/dx;

                denom_x1_J22w_x(i,j)=(Phi2(i,j)-Phi2(im,j))/dx;
                denom_x1_J22w_y(i,j)=((Phi2(i,jp)-Phi2(i,jm))/(2*dy)+(Phi2(im,jp)-Phi2(im,jm))/(2*dy))/2;

                denom_x1_J22(i,j)=sqrt(denom_x1_J22w_x(i,j)^2+denom_x1_J22w_y(i,j)^2);

                if(denom_x1_J22(i,j)==0)
                    denom_x1_J22(i,j)=tol;
                else
                    denom_x1_J22(i,j)=sqrt(denom_x1_J22w_x(i,j)^2+denom_x1_J22w_y(i,j)^2);
                end
                
                x1_J22(i,j)=x1_J22_1(i,j)*numer_x1_J22(i,j)/denom_x1_J22(i,j);

                %J22(i,j+1/2)
                
                y2_J22_1(i,j)=(J22_1(i,j)+J22_1(i,jp))/2;
                numer_y2_J22(i,j)=(Phi2(i,jp)-Phi2(i,j))/dy;

                denom_y2_J22n_x(i,j)=((Phi2(ip,j)-Phi2(im,j))/(2*dx)+(Phi2(ip,jp)-Phi2(im,jp))/(2*dx))/2;
                denom_y2_J22n_y(i,j)=(Phi2(i,jp)-Phi2(i,j))/dy;

                denom_y2_J22(i,j)=sqrt(denom_y2_J22n_x(i,j)^2+denom_y2_J22n_y(i,j)^2);

                if(denom_y2_J22(i,j)==0)
                    denom_y2_J22(i,j)=tol;
                else
                    denom_y2_J22(i,j)=sqrt(denom_y2_J22n_x(i,j)^2+denom_y2_J22n_y(i,j)^2);
                end
                
                y2_J22(i,j)=y2_J22_1(i,j)*numer_y2_J22(i,j)/denom_y2_J22(i,j);

                %J22(i,j-1/2);

                y1_J22_1(i,j)=(J22_1(i,j)+J22_1(i,jm))/2;
                numer_y1_J22(i,j)=(Phi2(i,j)-Phi2(i,jm))/dy;

                denom_y1_J22s_x(i,j)=((Phi2(ip,j)-Phi2(im,j))/(2*dx)+(Phi2(ip,jm)-Phi2(im,jm))/(2*dx))/2;
                denom_y1_J22s_y(i,j)=(Phi2(i,j)-Phi2(i,jm))/dy;

                denom_y1_J22(i,j)=sqrt(denom_y1_J22s_x(i,j)^2+denom_y1_J22s_y(i,j)^2);

                if(denom_y1_J22(i,j)==0)
                    denom_y1_J22(i,j)=tol;
                else
                    denom_y1_J22(i,j)=sqrt(denom_y1_J22s_x(i,j)^2+denom_y1_J22s_y(i,j)^2);
                end
                
                
                y1_J22(i,j)=y1_J22_1(i,j)*numer_y1_J22(i,j)/denom_y1_J22(i,j);







                % J22(i,j)=((pi*epsilon)/4)*(((hp_Phi2(i,j)*(1-hp_Phi2(i,j)))/sqrt(Phi2(i,j)*(1-Phi2(i,j)))))* Phi2_dot(i,j)./dt*(((g_Phi23(i,j))/((sqrt(g_Phi22(i,j)))*sqrt(g_Phi33(i,j)))))*...
                %             (((Com32(i,j)-Com22(i,j))*grad_x_Phi2(i,j)+(Com32(i,j)-Com22(i,j))*grad_y_Phi2(i,j)))*(1/sqrt(g_Phi22(i,j)));



             end


             

               

    end
end


for i=1:Nx 
    for j=1:Ny
         jp=j+1;
         jm=j-1;
         ip=i+1;
         im=i-1;
       
                %gradient J11

                grad_J11(i,j)=(x2_J11(i,j)-x1_J11(i,j))/dx+(y2_J11(i,j)-y1_J11(i,j))/dy;

               
                %gradient J12

                grad_J12(i,j)=(x2_J12(i,j)-x1_J12(i,j))/dx+(y2_J12(i,j)-y1_J12(i,j))/dy;


                
                %gradient J21

                grad_J21(i,j)=(x2_J21(i,j)-x1_J21(i,j))/dx+(y2_J21(i,j)-y1_J21(i,j))/dy;

                

                %gradient J22

                grad_J22(i,j)=(x2_J22(i,j)-x1_J22(i,j))/dx+(y2_J22(i,j)-y1_J22(i,j))/dy;









                



    end
end


     
     for i=1:Nx
         for j=1:Ny
            m(i,j)=a(i,j)*d(i,j)-b(i,j)*c(i,j);
            m1(i,j)=d(i,j)/m(i,j);
            m2(i,j)=-b(i,j)/m(i,j);
            m3(i,j)=-c(i,j)/m(i,j);
            m4(i,j)=a(i,j)/m(i,j);
            
            miuu1(i,j)=m1(i,j)*(fenzi1(i,j)-fenzi2(i,j)- grad_J11(i,j)-grad_J12(i,j))+m2(i,j)*(fenzi12(i,j)-fenzi22(i,j)-grad_J22(i,j)-grad_J21(i,j));
            miuu2(i,j)=m3(i,j)*(fenzi1(i,j)-fenzi2(i,j)- grad_J11(i,j)-grad_J12(i,j))+m4(i,j)*(fenzi12(i,j)-fenzi22(i,j)-grad_J22(i,j)-grad_J21(i,j));    
         end
     end
     
      miu(1:Nx,1:Ny)=miu(1:Nx,1:Ny)+(dt*miuu1(1:Nx,1:Ny));
      miu2(1:Nx,1:Ny)=miu2(1:Nx,1:Ny)+(dt*miuu2(1:Nx,1:Ny));

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Update Concentration Field%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Com11(1:Nx,1:Ny)=((miu(1:Nx,1:Ny)-O1+2*C1)/2/C1-(miu2(1:Nx,1:Ny)-P1+2*C1)/(B1+C1)/2)/((A1+C1)/C1-C1/(B1+C1));
        Com21(1:Nx,1:Ny)=((miu(1:Nx,1:Ny)-O2+2*C2)/2/C2-(miu2(1:Nx,1:Ny)-P2+2*C2)/(B2+C2)/2)/((A2+C2)/C2-C2/(B2+C2));
        Com31(1:Nx,1:Ny)=((miu(1:Nx,1:Ny)-O3+2*C3)/2/C3-(miu2(1:Nx,1:Ny)-P3+2*C3)/(B3+C3)/2)/((A3+C3)/C3-C3/(B3+C3));
        Com12(1:Nx,1:Ny)=((miu(1:Nx,1:Ny)-O1+2*C1)/2/(A1+C1)-(miu2(1:Nx,1:Ny)-P1+2*C1)/2/C1)/(C1/(A1+C1)-(B1+C1)/C1);
        Com22(1:Nx,1:Ny)=((miu(1:Nx,1:Ny)-O2+2*C2)/2/(A2+C2)-(miu2(1:Nx,1:Ny)-P2+2*C2)/2/C2)/(C2/(A2+C2)-(B2+C2)/C2);
        Com32(1:Nx,1:Ny)=((miu(1:Nx,1:Ny)-O3+2*C3)/2/(A3+C3)-(miu2(1:Nx,1:Ny)-P3+2*C3)/2/C3)/(C3/(A3+C3)-(B3+C3)/C3);
        Com13(1:Nx,1:Ny)=(1-Com11(1:Nx,1:Ny)-Com12(1:Nx,1:Ny));
        Com23(1:Nx,1:Ny)=(1-Com21(1:Nx,1:Ny)-Com22(1:Nx,1:Ny));
        Com33(1:Nx,1:Ny)=(1-Com31(1:Nx,1:Ny)-Com32(1:Nx,1:Ny));
      
      
      Com1(1:Nx,1:Ny)=Com11(1:Nx,1:Ny).*hp_Phi1(1:Nx,1:Ny)+Com21(1:Nx,1:Ny).*hp_Phi2(1:Nx,1:Ny)+Com31(1:Nx,1:Ny).*hp_Phi3(1:Nx,1:Ny);
      Com2(1:Nx,1:Ny)=Com12(1:Nx,1:Ny).*hp_Phi1(1:Nx,1:Ny)+Com22(1:Nx,1:Ny).*hp_Phi2(1:Nx,1:Ny)+Com32(1:Nx,1:Ny).*hp_Phi3(1:Nx,1:Ny);
      
 if((mod(istep,nprint)==0)|| (istep == 1))
        Phi0(1 : Nx, 1 : Ny, 1) = Phi1(1 : Nx, 1 : Ny);
        Phi0(1 : Nx, 1 : Ny, 2) = Phi2(1 : Nx, 1 : Ny);
        Phi0(1 : Nx, 1 : Ny, 3) = Phi3(1 : Nx, 1 : Ny);
   fprintf('done step: %5d\n',istep);
   toc;
   expr=['data_' num2str(istep) ' =istep; ' ];
   eval(expr)
   name=['data_' num2str(istep)];
   save(name,'Phi0' ,'Phi1','Phi2','Phi3','Phi','miu','ttime','miu2','Com1','Com2','Com11','Com12','Com21','Com22','Com31','Com32','hp_Phi1','hp_Phi2','hp_Phi3','Phi1_dot','Phi2_dot','Phi3_dot','Phi_dot','termg33','termh33','termi33','Term3','Term2','Term1','GG1','GG2','GG3','g1_3')

 end
end