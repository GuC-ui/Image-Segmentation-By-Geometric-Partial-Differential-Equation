clear

I = imread("Sample.png"); %画像を行列として読み込む。
%I = rgb_to_luminance(I);
%もし、画像がカラーならば、コメントアウトを消し、白黒画像に変換する。
I = 255 - I;
imagesc(I);colormap(gray(64)); %display the image 
%%%%%


syms a b c t n g(a,b,t) phi(t) 
g(a,b,t) = [a*cos(t) b*sin(t)]';

epsilon = 0.1;
%phi(t) = 1 - epsilon + epsilon*sqrt(1 - epsilon + epsilon*(t^2));
phi(t) = 1;
phi_(t) = diff(phi(t));

N = 50; %離散化を行う数
z = zeros(2,N+1); %x座標 : z(1,:), y座標 : z(2,:) 
lambda = 0.1;
time = 100; %経過時間
tau = 0.3; 
omega = 10;
[v,w] = size(I); %行列画像の行数と列数を、vとｗとして格納する。
Fmax = 45; %Fmax 
Fmin = -46;
G = 10;
by = 2.2;
a = 3.0;

m = 1;
r = 200; %初期曲線の半径
for i = 0+2*pi/N:2*pi/N:2*pi+2*pi/N
    z(:,m) = g(r,r,i)+[v/2;w/2]; %初期曲線の離散化
    m = m + 1;
end

for q = 1:time
   
    %The discretized arc-lngth 
    r = zeros(1,N+1);
    r(2:N) = ((z(1,2:N)-z(1,1:N-1)).^2 + (z(2,2:N)-z(2,1:N-1)).^2).^(1/2);
    r(1) = norm(z(:,1)-z(:,N));
    r(N+1) = r(1);
    
    %the average arc-lgth, which will be used for adjusted tangential
    %velocity
    rd = zeros(1,N);
    rd(1:N-1) = (r(1:N-1) + r(2:N))/2;
    rd(N) = (r(N) + r(1))/2;
    
    %tangential vector
    t = zeros(2,N+1);
    for i = 1:2
        t(i,2:N) = (z(i,2:N)-z(i,1:N-1))./r(2:N); 
    end
    t(:,1) = (z(:,1)-z(:,N))./r(1);
    t(:,N+1) = t(:,1);
    
    %tangential angle
    nu = zeros(1,N+3);
    
    if t(2,1) < 0
        nu(2) = 2*pi - acos(t(1,1));
    else
        nu(2) = acos(t(1,1));
    end
    
    for i = 2:N+1
        if dot(t(:,i-1),t(:,i)) > 0
            nu(i+1) = nu(i) + acos(det([t(:,i-1)';t(:,i)']));  
        elseif det([t(:,i-1)';t(:,i)']) > 0
            nu(i+1) = nu(i) + acos(dot(t(:,i-1),t(:,i)));
        else
            nu(i+1) = nu(i) - acos(dot(t(:,i-1),t(:,i)));
        end
    end
    nu(1) = nu(2) - (nu(N+2) - nu(N+1));
    nu(N+3) = nu(N+2) + (nu(3) - nu(2));
    
    % average normal vector
    avt = zeros(2,N); %averaged tangential vector 
    nd = zeros(2,N); 
    for i = 1:2
        avt(i,1:N-1) = (t(i,1:N-1) + t(i,2:N))./2;
    end
    avt(:,N) = (t(:,1) + t(:,N))./2;
    avt(:,1:N) = avt(:,1:N)./norm(avt(:,1:N));
    nd(:,1:N) = [0 -1;1 0]*avt(:,1:N); %rotated by the apropriate matrix
  
    %Curvature
    k = zeros(1,N);
    k(1:N) = (nu(3:N+2)-nu(1:N))./(2*r(1:N));
    
    kd = zeros(1,N);  
    kd(1:N-1) = (k(2:N) + k(1:N-1))/2;
    kd(N) = (k(1) + k(N))/2;
    
    %External Force
    F = zeros(1,N+1); 
    xtmp = ceil(z(1,1:N));
    ytmp = ceil(z(2,1:N));
    for i = 1:N
        F(i) = Fmax - (Fmax - Fmin)*(sigmoid(double(I(xtmp(i),ytmp(i)))/255,a)-sigmoid(0,a))/(sigmoid(1,a)-sigmoid(0,a));
    end
    F(N+1) =F(1);
   
    
    %phi = 1,omega = 0 : asymptotic uniform distribution method 
    ppbeta = zeros(1,N); %the second derivative of beta
    tic
    ppbeta(2:N-1) = (kd(3:N)+F(3:N)-(kd(2:N-1)+F(2:N-1)))./(0.5*r(2:N-1).*(r(2:N-1)+r(3:N))) - (kd(2:N-1)+F(2:N-1)-(kd(1:N-2)+F(1:N-2)))./(0.5*r(2:N-1).*(r(1:N-2)+r(2:N-1)));
    ppbeta(1) = (kd(2)+F(2)-(kd(1)+F(1)))/(0.5*r(1)*(r(1)+r(2))) - (kd(1)+F(1)-(kd(N)+F(N)))/(0.5*r(1)*(r(N)+r(1)));
    ppbeta(N) = (kd(1)+F(1)-(kd(N)+F(N)))/(0.5*r(N)*(r(N)+r(1))) - (kd(N)+F(N)-(kd(N-1)+F(N-1)))/(0.5*r(N)*(r(N-1)+r(N)));
    
    f = (ppbeta + (k.^2).*(kd+F(1:N))).*phi_(k)-k.*(kd+F(1:N)).*phi(k);
    L = sum(r(1:N)); %total length of the curve
    avef = dot(f,r(1:N))/L; %<f>:the average f
    avephi = dot(phi(k),r(1:N))/L;%<phi>:the average phi 
    psi = (avef.*r(1:N).*phi(k)./avephi) - (f.*r(1:N))+(L*avephi/N - phi(k).*r(1:N))*omega;
   
    Psi = zeros(1,N);
    for i = 2:N
        Psi(i) = Psi(i-1)+psi(i);
    end
    %Psi(2:N) = Psi(1:N-1) + psi(2:N);
    
    %Adjusted tangential velocity
    alpha = zeros(1,N);
    alpha(1) = alpha(1) - (dot(Psi,rd)-Psi(1)*rd(1))/(L * phi(kd(1))); 
    alpha(2:N) = (phi(kd(1))*alpha(1) + Psi(2:N))./phi(kd(2:N));
    %%%Semi-implicit numerical scheme%%%
    a = alpha./(2*rd);
    b = 1./rd;
    
    aplus = zeros(1,N);
    aplus(1:N) = b(1:N)/r(2:N+1) + a(1:N);
    aminus = b./r(1:N) - a;
    a0 = -(aminus + aplus);
    abal = abs(alpha); %the absolute value of the velocity
    
    tau = min(r)/((4 + 4*lambda)*(1/min(r) + max(abal)/2)); %Time step
   
    X = zeros(N,N); %semi-impliit numerical schemeのための変換行列
    for i = 2:N-1
          X(i-1,i) = -aminus(i)*tau;
          X(i,i) = (1 - a0(i)*tau);
          X(i+1,i) = -aplus(i)*tau;
    end
    X(1,1) = 1 - a0(1)*tau;
    X(2,1) = -aplus(1)*tau;
    X(N,1) = -aminus(1)*tau;
    X(1,N) = -aplus(N)*tau;
    X(N-1,N) = -aminus(N)*tau;
    X(N,N) = 1 - a0(N)*tau;

    Ztmp = zeros(2,N);
    for i = 1:2
        Ztmp(i,1:N) = z(i,1:N) + nd(i,1:N).*tau.*F(1:N); 
    end   
    z = Ztmp/X;%時間変化させるため、連立一次方程式を解く。
    %%%%%%%%%%%%%%%%%%%%%%
    
    z(:,N+1) = z(:,1); %周期条件として、曲線の端点を一致させる。
 
    z = [0 1; -1 0]*z;
    z = [1 0; 0 -1]*z;
 
   if q == time
       I = 255 - I;
       figure
       hold on
       imagesc(I);colormap(gray(64));
       plot(z(1,:),z(2,:),'LineWidth',1.5)
   end

   z = [1 0;0 -1]*z;
   z = [0 -1;1 0]*z;
end

 
