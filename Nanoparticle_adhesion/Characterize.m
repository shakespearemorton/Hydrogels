% create a 2D grid
clc
clear all
Atom=zeros(1,4);
row=1;
l = REPLACEL;
m = REPLACEM;
s = REPLACES;
rfactor = REPLACER;
th = linspace(0,pi,200);    % inclination
phi = linspace(0,2*pi,200); % azimuth
[th,phi] = meshgrid(th,phi);
Y = s*(rfactor+harmonicY(l,m,th,phi,'type', 'real'));
r = abs(Y);
maxi = max(max(r));
r = r /(maxi/10);
s = 1/(maxi/10);
Y = s*(rfactor+harmonicY(l,m,th,phi,'type', 'real'));
r = abs(Y);
maxi = max(max(r));
save=r;
hei = r;
[x,y,z] = sph2cart(phi,pi/2-th,r);
%colormap(winter)
%shading interp
surf(x,y,z)
[w,n] = size(z);
area1 = 0;
for i = 1:w-1
      for j = 1:n-1
          v0 = [x(i,j)     y(i,j)     z(i,j)    ];
          v1 = [x(i,j+1)   y(i,j+1)   z(i,j+1)  ];
          v2 = [x(i+1,j)   y(i+1,j)   z(i+1,j)  ];
          v3 = [x(i+1,j+1) y(i+1,j+1) z(i+1,j+1)];
          a = v1 - v0;
          b = v2 - v0;
          c = v3 - v0;
          A = 1/2*(norm(cross(a, c)) + norm(cross(b, c)));
          area1 = area1 + A;
      end
end
x1 = x(:);
y1 = y(:);
z1 = z(:);
P=[x1 y1 z1];
%shp = alphaShape(P);
%SA = surfaceArea(shp);
%plot(shp)
Atom(row,1)=area1;
[K,H,Pmax,Pmin,K2,H2] = surfature(x,y,z);
d1=3.14/199;
d2=6.28/199;
x = 0:d1:3.14;
y = 0:d2:6.28;
Ih = trapz(y,trapz(x,H2,2))/(area1);
Atom(row,2)=Ih;
r = r-min(min(r));
zmax = max(max(r));
rms = sqrt(sum(sum(r.^2))/(200*200));
Atom(row,3)=rms;
hei = hei(:);
Ra = sum(hei)/(200*200);
wavy = 0;
for j = 2:(200*200)
    wavy = wavy+(abs((hei(j)-hei(j-1)))/(pi/100));
end
wavy = wavy/((200*200)-1);
avg_wave = (2*pi*Ra)/wavy;
Atom(row,4)=avg_wave;
dlmwrite('_Characterize.txt', Atom)


