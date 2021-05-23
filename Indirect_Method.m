%3次元の曲面（z=x^2+y^2）を偏微分方程式により時間変化させ、その等高線を
%画像の輪郭に変える事ができる。

clear

I = imread("pokemon.png"); 
%画像を、行列として読み込む。
%なお、画像がカラーならば、行列の要素毎に0～255の数が３つ入っており、（6次元）
%画像が白黒ならば、0～255の数が１つ入っている。（2次元）

I_ = I;

I = rgb2lab(I); %もし画像がカラーならば、分析をするために
I = I(:,:,1);   %白黒の画像に変換する

[v,w] = size(I); %行列の行数、列数を格納

x = 0:1:v-1;
y = 0:1:w-1; %離散化する

[X,Y] = meshgrid(x,y);
X = transpose(X);
Y = transpose(Y);
r = 1;
h1 = 1; h2 = 1;
time = 100; %どれだけ時間を変化させるか。
F = -2*((X-v/2).^2 + (Y-w/2).^2) + (min(v/2-5,w/2-5))^2;
%画像行列に合わせて、z=x^2+y^2を変える。
delta = 2^(-50);
sigma = 2;
tau = 0.1;

a = 8; %シグモイド関数のパラメータ－
Fmax = 30; Fmin = -100; %外力のパラメーター

ve = [0,0];

%%フィルター処理を行うためのパラメーター
KD1F = [0 -1 0;0 0 0;0 1 0];
KD2F = [0 0 0; -1 0 1;0 0 0];

KDp1F = [0 0 0;0 -1 0;0 1 0];
KDp2F = [0 0 0;0 -1 1; 0 0 0];
KDm1F = [0 1 0;0 -1 0;0 0 0];
KDm2F = [0 0 0;1 -1 0; 0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:time
  %%フィルター処理を行う。
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
  
  %%%%%%%偏微分方程式であらわされる外力による時間変化%%%%%%%%
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
