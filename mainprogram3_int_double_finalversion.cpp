//
//  mainprogram3_int_double.cpp
//  testpr3paper
//


#include<stdio.h>
#include <math.h>
#include "int_double12.0.h"

long floor(int_double a)
{
  long l,r;
  l=a.left;
  r=-a.right;
  if(l==r)
    return l;
  print_int_double_str("Fatal error in floor with ",a);
  exit(0);
  return 0;
}

long cheat_floor(int_double a)
{
  return -a.right;
}

bool ge(int_double a, int_double b)
{
  if(a.left>=-b.right)
    return true;
  if(-a.right<b.left)
    return false;
  printf("Overlap in ge. Exiting.\n");
  exit(0);
  return 0;
}

#define le(a,b) ge(b,a)

bool gt(int_double a, int_double b)
{
 
  if(a.left>-b.right)
    return true;
  if(-a.right<=b.left)
    return false;
  printf("Overlap in gt. Exiting.\n");
  exit(0);
  return 0;
}
#define lt(a,b) gt(b,a)



int_double iunion(int_double a, int_double b)
{
  int_double res;
  if(a.left<=b.left)
    res.left=a.left;
  else
    res.left=b.left;
  if(-a.right>=-b.right)
    res.right=a.right;
  else
    res.right=b.right;
  return res;
}

long k,g,h,s,r,t,g0,h0,g1,h1,flag;
int_double mu1,mu2,xi,lam,lam1,lam2,D,sigma,Y,goal;
void calc(int_double *ex, int_double *c, int pr)
{

    int_double kk,logk,k2;
    int_double th,rr,ss,tt,gg,hh,rho,H,E1,E2,E3,m1,m2,Z0,Z1,
           reta,logC1,logC2,logC3,logC;
    k=floor(lam/(1.0-mu1-mu2)+0.000003);
    if(k<129)
      {
	printf("k too small. Exiting.\n");
	exit(0);
      }    kk=(double) k;
    logk=log(kk); k2=k*k;
    rho=int_double(320862,320864); th=int_double(217719,217721); rho/=10000;th/=10000;
    if (k<=89999) {rho=int_double(3.205500,3.205501); th=int_double(1.77774,1.77775);}
    if (k<=499) {rho=int_double(3.1964972028,3.1964972029); th=int_double(2.2435101628,2.2435101629);}
    if (k<=339) {rho=int_double(3.1929499395,3.1929499396); th=int_double(2.3331208028,2.3331208029);}
    if (k<=190) {rho=int_double(3.1841274238,3.1841274239); th=int_double(2.3533307441,2.3533307442);}
    if (k<=170) {rho=int_double(3.1818685121,3.1818685122); th=int_double(2.3792823257,2.3792823258);}
    if (k<=148) {rho=int_double(3.1788714390,3.1788714391); th=int_double(2.3825812477,2.3825812478);}
    if (k<=146) {rho=int_double(3.1785513229,3.1785513230); th=int_double(2.3916650550,2.3916650551);}
    if (k<=139) {rho=int_double(3.1775270431,3.1775270432); th=int_double(2.3952838649,2.3952838650);}
    if (k<=137) {rho=int_double(3.1772070968,3.1772070969); th=int_double(2.4092926606,2.4092926607);}

    
    r=cheat_floor(rho*k2+1.0);
   
    rr=r; ss=s;
    gg=g; hh=h; tt=t;
    m1=floor(lam/(1.0-mu1));
    m2=floor(lam/(1.0-mu2));
    Z0=0.5*((m1*m1+m1)*(1.0-mu1)+(m2*m2+m2)*(1.0-mu2)-
        -hh*hh+hh-(1.0-mu1-mu2)*(gg*gg+gg));
    Z1=hh+gg-m1-m2-1.0;
    if (Z1.left<0.0)
      {
	if(-Z1.right>=0.0)
	  {
	    print_int_double_str("Z1 indeterminate in calc. Exiting.",Z1);
	    exit(0);
	  }
	else
	  H=Z0+lam2*Z1;
      }
    else
      H=Z0+lam1*Z1;
    reta=xi*pow(gg,1.5);
    E1=0.001*k2;
    E2=0.5*tt*(tt-1.0)+hh*tt*exp(-ss/(hh*tt))+ss*ss/(2.0*tt*reta);
    E3=log(Y*lam1*lam1)/(1.0*Y*lam1*lam1*lam1*lam1);
    *ex=(-E3+(1.0/(2.0*rr*ss))*(H-mu1*E1-mu2*E2))*lam1*lam1;
    logC1=th*k2*kk*logk;
    logC2=ss*ss/tt+10.5*xi*xi*tt*gg*gg*log(gg)*log(gg)/D;
    logC2 -=(ss*log(0.1*reta)*((reta+hh)*pow(1.0-1.0/hh,ss/tt)-h));
    logC3=1.0417*reta*log(10.4167*reta);
    logC=logC3/rr+(4.65*lam2*log(lam2)+logC1+logC2)/(2.0*rr*ss);
    *c=exp(logC)+1.0/kk;
    if (pr==1){
      printf("%8.4f-%8.4f %4ld",lam1.left,-lam2.right,k);
        if(g>0) {printf("%3ld %2ld %2ld %2ld ",
			s,g-g0,h1-h,t);
	  print_int_double(1.0/(*ex)+0.00005);
	  print_int_double_str(" ",*c+0.00005);
	}
        else printf("\n");
    }
   
}

int main(int argc, char**argv)
{
  _fpu_rndd();
  printf("Command:-");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=8)
    {
      printf("Usage:- %s <Y> (<xi> <xi>) (<sigma> <sigma>) lam8 lam9\n",argv[0]);
      return 0;
    }

  Y=atol(argv[1]);
  xi=int_double(atof(argv[2]),atof(argv[3]));
  sigma=int_double(atof(argv[4]),atof(argv[5]));
  long lam8,lam9;
  lam8=atol(argv[6]);lam9=atol(argv[7]);
  if ((lam9<lam8)||(lam8<=80)||(lam9>=300))
    {
      printf("need lam8<=lam9 and lam8>=80 and lam9<300. Exiting.\n");
      return 0;
    }
  
    int_double E,r[9],tmp,maxex,con,maxcon;
    int_double bestth,bestcon,bp[5000];
    long i,j,i0,w,n,m,bestg,besth,bests,s0,s1;
    mu1=1905;mu1/=10000; mu2=1603;mu2/=10000;
    goal=13294357;goal/=100000;
    D=1019*Y;D/=10000;
    if (sigma.left<0.0) flag=1; else flag=0;

    bp[1]=lam8; bp[2]=lam9; j=3;
    i0=floor(lam9/(1.0-mu1-mu2))+10;
    for (i=1; i<=i0;i++)
      {
	w=i;
	r[1]=w*(1.0-mu1);
	r[2]=w*(1.0-mu2);
	r[3]=(w-0.000003)*(1.0-mu1-mu2);
	for (m=1;m<=3;m++)
	  if (lt(r[m],lam9) && gt(r[m],lam8))
	    bp[j++]=r[m];
      }
    n=j-1;
    for (i=1;i<=n-1;i++)
      for(j=i+1;j<=n;j++)
	if (lt(bp[j],bp[i]))
	  {tmp=bp[i]; bp[i]=bp[j]; bp[j]=tmp;}
    maxex=0.0;
    maxcon=0.0;
    for (j=1;j<=n-1;j++)
      {
	
	lam=0.5*(bp[j]+bp[j+1]);
	lam1=bp[j]; lam2=bp[j+1];
	g0=floor(lam/(1.0-mu1)+1.0); g1=g0+1;
	h1=floor(lam/(1.0-mu2)); h0=h1-1;
	bestg=-1; besth=-1; bestth=1.0e20; bestcon=1.0e40;
	for (g=g0;g<=g1;g++)
	  for (h=h0;h<=h1;h++)
	    {
	      t=g-h+1;
	    
	      if ((g>=100) && (le(g*1000,1254*lam1)))
		{
		  if (flag==0)
		    {
		      s0=floor(sigma*h*t+1.0);
		      s1=s0;
		    }
		  else
		    {
		      s0=h*(t-1)/4;
		      s1=h*t/2;
		    }
		  for(s=s0; s<=s1; s++) {
		    calc(&E,&con,0);
		    int_double tmp=1.0/E;
		    if((E.left>0.0)&&(-tmp.right<goal.left)&&(-con.right<bestcon.left))
		      {
			bestth=tmp; bestg=g; besth=h; bests=s;
		      }
		  }
		}
	    }
        
      g=bestg; h=besth; t=g-h+1;
      s=bests;
      calc(&E,&con,0);
      int_double tmp=1.0/E;
      maxex=iunion(maxex,tmp);
      maxcon=iunion(maxcon,con); 
      printf("lam = %f %f: ex: %10.8f const: %10.8f\n",lam1.left,-lam2.right,-tmp.right,-con.right);
    }
    printf("Final ex       = %14.12e\n",-maxex.right);
    printf("Final constant = %14.12e\n",-maxcon.right);
    return 0;
}
