function [C_aug,SU_zeros]=Find_Square_Up_v1(A,B,C,des_zeros)
%%%%%%% In courtesy of Max Zheng Qu, 05/02/2013%%%%%%%%%%%
%%%%%%% This function implement Misra's method on square up C matrix
%%%%%%% The code works for find a augmented C_aug to place the system zeros
%%%%%%% at des_zeros
%%%%%%% To find a augmented B_aug, one needs to transpose everything and
%%%%%%% treat B=C' and C=B'. The dual identity of the system  %%%%%%%%%%%%


OL_name=['{',inputname(1),',',inputname(2),',',inputname(3),'}'];
num_state=length(A);
num_input=size(B,2);
num_output=size(C,1);

D=zeros(num_output,num_input);
OL_tzeros=tzero(A,B,C,D);

if length(OL_tzeros)==0;
disp(strcat(['There are no tzeros in the original system']));
else
disp(strcat(['There are tzeros in the original system'],[' ','they are',': '],num2str(OL_tzeros)));
end

% C_bar=[C;1,0,0,0,0,0];
% D_bar=[D;0,0,0];
% SU_tzeros=tzero(A,B,C_bar,D_bar);


B_perp=(null(B'));

T=[B';B_perp'];

At=T*A*inv(T);
Bt=T*B;
Ct=C*inv(T);
Dt=D;

Bt_0=sum(Bt,2);

index_0=find(abs(Bt_0)>1e-8);
n1=length(index_0);
n2=num_state-n1;


A22=At(n1+1:n1+n2,n1+1:n1+n2);
A21=At(n1+1:n1+n2,1:n1);
C11=Ct(:,1:n1);
C12=Ct(:,n1+1:n1+n2);

C21=(null(C11))';
C1=[C11;C21];

C2tilt=[C12;zeros(num_input-num_output,num_state-num_input)];

A22tilt=A22-A21*inv(C1)*C2tilt;
% des_zeros=[-0.1,-0.2,-0.3];

B_pseudo=A21*inv(C1);
B_pseudo_2=B_pseudo(:,n1-num_input+num_output+1:n1);

num_zeros=length(A22tilt);
disp(strcat(['Square up can place:',' '],num2str(num_zeros),[' ', 'zeros']));
C2hat_p=place(A22tilt,B_pseudo_2,des_zeros);
C2hat=[zeros(num_output,num_state-num_input);C2hat_p];
C2=C2tilt+C2hat;

A22tiltFB=A22tilt-A21*inv(C1)*C2hat;
% eig_A22tiltFB=Check_UnstaPole_v1(A22tiltFB);
A22FB=A22-A21*inv(C1)*C2;
% eig_A22FB=Check_UnstaPole_v1(A22FB);

Ct_aug=[C1,C2];

C_aug=Ct_aug*T;

D_aug=[D;zeros(num_input-num_output,num_input)];

SU_zeros=tzero(A,B,C_aug,D_aug);


disp(strcat(['There are tzeros in the augmented system'],[' ','they are',': '],num2str(SU_zeros)));





end