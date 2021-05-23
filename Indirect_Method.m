%3�����̋Ȗʁiz=x^2+y^2�j��Δ����������ɂ�莞�ԕω������A���̓�������
%�摜�̗֊s�ɕς��鎖���ł���B

clear

I = imread("pokemon.png"); 
%�摜���A�s��Ƃ��ēǂݍ��ށB
%�Ȃ��A�摜���J���[�Ȃ�΁A�s��̗v�f����0�`255�̐����R�����Ă���A�i6�����j
%�摜�������Ȃ�΁A0�`255�̐����P�����Ă���B�i2�����j

I_ = I;

I = rgb2lab(I); %�����摜���J���[�Ȃ�΁A���͂����邽�߂�
I = I(:,:,1);   %�����̉摜�ɕϊ�����

[v,w] = size(I); %�s��̍s���A�񐔂��i�[

x = 0:1:v-1;
y = 0:1:w-1; %���U������

[X,Y] = meshgrid(x,y);
X = transpose(X);
Y = transpose(Y);
r = 1;
h1 = 1; h2 = 1;
time = 100; %�ǂꂾ�����Ԃ�ω������邩�B
F = -2*((X-v/2).^2 + (Y-w/2).^2) + (min(v/2-5,w/2-5))^2;
%�摜�s��ɍ��킹�āAz=x^2+y^2��ς���B
delta = 2^(-50);
sigma = 2;
tau = 0.1;

a = 8; %�V�O���C�h�֐��̃p�����[�^�|
Fmax = 30; Fmin = -100; %�O�͂̃p�����[�^�[

ve = [0,0];

%%�t�B���^�[�������s�����߂̃p�����[�^�[
KD1F = [0 -1 0;0 0 0;0 1 0];
KD2F = [0 0 0; -1 0 1;0 0 0];

KDp1F = [0 0 0;0 -1 0;0 1 0];
KDp2F = [0 0 0;0 -1 1; 0 0 0];
KDm1F = [0 1 0;0 -1 0;0 0 0];
KDm2F = [0 0 0;1 -1 0; 0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:time
  %%�t�B���^�[�������s���B
  Dp1F = imfilter(F,KDp1F);
  Dp2F = imfilter(F,KDp2F);
  Dm1F = imfilter(F,KDm1F);
  Dm2F = imfilter(F,KDm2F);
  D1F = imfilter(F,KD1F);
  D2F = imfilter(F,KD2F);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  G = sqrt((abs(Dp1F)+abs(Dm1F)).^2 + (abs(Dp2F)+abs(Dm2F)).^2)./2;
  div = imfilter(D1F./((G.^sigma)+delta).^(1/sigma),KD1F) + imfilter(D2F./((G.^sigma)+delta).^(1/sigma),KD2F);
  
  %E = zeros(size(I));
  
  %%%%%%%�Δ����������ł���킳���O�͂ɂ�鎞�ԕω�%%%%%%%%
  F = F + tau*G.*(div - (Fmax - (Fmax - Fmin)*(sigmoid(double(I)/255,a)-sigmoid(0,a))/(sigmoid(1,a)-sigmoid(0,a))));
  %%%%%%%%%%%%%%%%%%%
  %F = F + tau* G.*max(Fmax - (Fmax - Fmin)* double(I)/255,Fmax - (Fmax - Fmin)*(sigmoid(double(I)/255,a)-sigmoid(0,a))/(sigmoid(1,a)-sigmoid(0,a)));
  DF = imfilter(F,KD2F).^2 + imfilter(F,KD1F).^2;
  
  if(mod(k,3) == 0)
  F = F./sqrt(F.^2 + DF*1);
  end
end

I_ = rot90(I_,-1);
I_ = fliplr(I_);
figure
hold on
imagesc(I_);colormap(gray(64));
contour(X,Y,F,ve,'r','LineWidth',1.0);
figure 
surface(X,Y,F);
