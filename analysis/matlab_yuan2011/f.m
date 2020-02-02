%Data structure
%RIXIndicator trait  parent1 parent2 RIXmarker1-markerP
function hatab=fsample(data)
aaa = clock;
drop=1000;
%drop the number for dropping out in the posterior sample
delta=0.001;
%delta prior parameter
thin=50;
%thin the thinning number for the posterior sample
hmax=8000;
%hmax the number of posterior sample
num.marker=length(data(1,:))-4;
n=length(data(:,1));
data34=data(:,3:4);
num.RI=max(data34(:));
%transform marker matrix x and w
x=data(:,5:end);
w=data(:,5:end);
x(x==1)=sqrt(2);
x(x==-1)=-sqrt(2);
w(w==1)=-1;
w(w==0)=1;
xsquare=sum(x.^2);
wsquare=sum(w.^2);
%A the matrix in the model for random effect
%A is from the parent information p1p2
A=zeros(n, num.RI);
for i=1:n
   A(i,data(i,3))=1;
   A(i,data(i,4))=1;
end
Asquare=sum(A.^2);
%trait value y
y=data(:,2);
resultu=zeros(1,hmax);
resultsigmaA=resultu;
resultsigma0=resultu;
resulta=zeros(hmax,num.marker);
resultb=resulta;
bold1=zeros(1,num.marker);
newa=bold1;
bold2=bold1;
newb=bold1;
newsigmasquare=bold1;
newvsquare=bold1;
bold3=1;
bold4=ones(1,num.marker);
bold5=ones(1,num.marker);
bold6=zeros(1,num.RI);
newalpha=bold6;
bold7=1;
z=[x,w,A];
%do the sampling
for h=1:hmax
boldu=[bold1,bold2,bold6];
meanu=mean(y-z*boldu');
newu=(sqrt(bold7/n))*randn+meanu;
bold1(1)=0;
bolda=[bold1,bold2,bold6];
meana=(xsquare(1)*bold4(1)+bold7)^(-1)*bold4(1)*sum(x(:,1).*(y-newu-z*bolda'));
vara=(xsquare(1)*bold4(1)+bold7)^(-1)*(bold7*bold4(1));
newa(1)=sqrt(vara)*randn+meana;
newsigmasquare(1)=newa(1)^2/chi2rnd(1-2*delta);
bold1(1)=newa(1);
bold2(1)=0;
boldb=[bold1,bold2,bold6];
meanb=(wsquare(1)*bold5(1)+bold7)^(-1)*bold5(1)*sum(w(:,1).*(y-newu-z*boldb'));
varb=(wsquare(1)*bold5(1)+bold7)^(-1)*bold7*bold5(1);
newb(1)=sqrt(varb)*randn+meanb;
newvsquare(1)=newb(1)^2/chi2rnd(1-2*delta);
for j=2:num.marker
bold1(1:(j-1))=newa(1:(j-1));
bold2(1:(j-1))=newb(1:(j-1));
bold1(j)=0;
bolda=[bold1,bold2,bold6];
meana=(xsquare(j)*bold4(j)+bold7)^(-1)*bold4(j)*sum(x(:,j).*(y-newu-z*bolda'));
vara=(xsquare(j)*bold4(j)+bold7)^(-1)*(bold7*bold4(j));
newa(j)=sqrt(vara)*randn+meana;
newsigmasquare(j)=newa(j)^2/chi2rnd(1-2*delta);
bold1(j)=newa(j);
bold2(j)=0;
boldb=[bold1,bold2,bold6];
meanb=(wsquare(j)*bold5(j)+bold7)^(-1)*bold5(j)*sum(w(:,j).*(y-newu-z*boldb'));
varb=(wsquare(j)*bold5(j)+bold7)^(-1)*bold7*bold5(j);
newb(j)=sqrt(varb)*randn+meanb;
newvsquare(j)=newb(j)^2/chi2rnd(1-2*delta);
end
boldsigma0=[newa,newb,bold6];
newsigma0square=sum((y-newu-z*boldsigma0').^2)/chi2rnd(n);
for k=1:num.RI
bold6(1:k)=newalpha(1:k);
bold6(k)=0;
boldalpha=[newa,newb,bold6];
meanal=(Asquare(k)*bold3+newsigma0square)^(-1)*bold3*sum(A(:,k).*(y-newu-z*boldalpha'));
varal=(Asquare(k)*bold3+newsigma0square)^(-1)*newsigma0square*bold3;
newalpha(k)=sqrt(varal)*randn+meanal;
end
newsigmaAsquare=(sum(newalpha.^2))/chi2rnd((num.RI-2*delta));
resultu(h)=newu;
resulta(h,:)=newa;
resultb(h,:)=newb;
resultsigma0(h)=newsigma0square;
resultsigmaA(h)=newsigmaAsquare;
bold1=newa;
bold2=newb;
bold3=newsigmaAsquare;
bold4=newsigmasquare;
bold5=newvsquare;
bold6=newalpha;
bold7=newsigma0square;
end
%drop out the drop number
ra=resulta((drop+1):end,:);
rb=resultb((drop+1):end,:);
ru=resultu((drop+1):end);
rsigma0=resultsigma0((drop+1):end);
rsigmaA=resultsigmaA((drop+1):end);
%get the posterior sample
pos=(1:((hmax-drop)/thin));
rra=ra(thin*pos,:);
rrb=rb(thin*pos,:);
rru=ru(thin*pos);
rrsigma0=rsigma0(thin*pos);
rrsigmaA=rsigmaA(thin*pos);
%get the posterior mean
hatu=mean(rru);
hatsigma0=mean(rrsigma0);
hatsigmaA=mean(rrsigmaA);
hata=mean(rra);
hatb=mean(rrb);
hatab=[hata;hatb];
axis=(0:(num.marker-1))*5;
%axis is genome location based on data 
subplot(2,1,1);
stem(axis,hata,'.-')
title('additive effect')
xlabel('genome location')
subplot(2,1,2);
stem(axis,hatb,'.-')
title('dominance effect')
xlabel('genome location')

bbb=clock;
% function to obtain the posterior sample and get the posterior mean
% hatab is a 2*num.marker matrix, the first row is additive effect, the
% second row is the dominance effect



