//
//  main.c
//  testpr3paper
//
//  Created by chiara bellotti on 16/1/2024.
//

#include<stdio.h>
#include <math.h>
long k,g,h,s,r,t,g0,h0,g1,h1,flag;
double mu1,mu2,xi,lam,lam1,lam2,D,sigma,Y,goal;
void calc(ex,c,pr)
    double *ex,*c; int pr;
{
    double kk,logk,k2,log(),exp(),floor(),ceil();
    double th,rr,ss,tt,gg,hh,rho,H,E1,E2,E3,m1,m2,Z0,Z1,
           reta,logC1,logC2,logC3,logC;
    k=(long) (lam/(1.0-mu1-mu2)+0.000003);
    kk=(double) k;
    logk=log(kk); k2=kk*kk;
    rho=3.20863; th=2.17720;
    if (k<=89999) {rho=3.205502; th=1.77775;}
    if (k<=499) {rho=3.196497; th=2.24352;}
    if (k<=339) {rho=3.192950; th=2.33313;}
    if (k<=190) {rho=3.184127; th=2.35334;}
    if (k<=170) {rho=3.181869; th=2.37929;}
    if (k<=148) {rho=3.178871; th=2.38259;}
    if (k<=146) {rho=3.178551; th=2.39167;}
    if (k<=139) {rho=3.177527; th=2.39529;}
    if (k<=137) {rho=3.177207; th=2.40930;}
    r=(long) (rho*k2+1.0);
    rr=(double) r; ss=(double) s;
    gg=(double) g; hh=(double) h; tt=(double) t;
    m1=floor(lam/(1.0-mu1));
    m2=floor(lam/(1.0-mu2));
    Z0=0.5*((m1*m1+m1)*(1.0-mu1)+(m2*m2+m2)*(1.0-mu2)-
        -hh*hh+hh-(1.0-mu1-mu2)*(gg*gg+gg));
    Z1=hh+gg-m1-m2-1.0;
    if (Z1<0.0) H=Z0+lam2*Z1;
    else H=Z0+lam1*Z1;
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
        printf("%8.4f-%8.4f %4ld",lam1,lam2,k);
        if(g>0) printf("%3ld %2ld %2ld %2ld %9.4f %7.4f\n",
                       s,g-g0,h1-h,t,1.0/(*ex)+0.00005,*c+0.00005);
        else printf("\n");
    }
}
int main()
{
    double E,lam8,lam9,r[9],tmp,maxex,con,maxcon;
    double bestth,bestcon,bp[5000];
    long i,j,i0,w,n,m,bestg,besth,bests,s0,s1;
    mu1=0.1905; mu2=0.1603;
    goal=132.94357;
    while (1) {
        printf("enter Y:"); scanf("%lf", &Y);
        D=0.1019*Y;
        printf("enter xi: "); scanf("%lf", &xi);
        printf("enter sigma :"); scanf("%lf", &sigma);
        if (sigma<0.0) flag=1; else flag=0;
        printf("enter lambda range: ");
        scanf("%lf %lf", &lam8, &lam9);
        if ((lam9<lam8)||(lam8<=80.0)||(lam9>=300.0)) continue;
        printf(" approx. \n");
        printf("lambda range    k   s   a   b   t   exp const\n");
        printf("-------------- --- --- --- --- --- -----------\n");
        bp[1]=lam8; bp[2]=lam9; j=3;
        i0=(long) (lam9/(1.0-mu1-mu2))+10;
        for (i=1; i<=i0;i++){
            w=(double) i;
            r[1]=w*(1.0-mu1);
            r[2]=w*(1.0-mu2);
            r[3]=(w-0.000003)*(1.0-mu1-mu2);
            for (m=1;m<=3;m++) if ((r[m]<lam9) && (r[m]>lam8))
            bp[j++]=r[m];
        }
        n=j-1;
        for (i=1;i<=n-1;i++) for(j=i+1;j<=n;j++)
        if (bp[j]<bp[i]) {tmp=bp[i]; bp[i]=bp[j]; bp[j]=tmp;}
        maxex=0.0;
        maxcon=0.0;
        for (j=1;j<=n-1;j++){
            lam=0.5*(bp[j]+bp[j+1]);
            lam1=bp[j]; lam2=bp[j+1];
            g0=(long) (lam/(1.0-mu1)+1.0); g1=g0+1;
            h1=(long) (lam/(1.0-mu2)); h0=h1-1;
            bestg=-1; besth=-1; bestth=1.0e20; bestcon=1.0e40;
            for (g=g0;g<=g1;g++) for (h=h0;h<=h1;h++){
                t=g-h+1;
                if ((g>=100) && ((double) g<=1.254*lam1)){
                    if (flag==0) {
                        s0=(long) (sigma*h*t+1.0); s1=s0;
                    }
                    else {
                        s0=h*(t-1)/4;
                        s1=h*t/2;
                    }
                    for(s=s0; s<=s1; s++) {
                        calc(&E,&con,0);
                        if((E>0.0)&&(1.0/E<goal)&&(con<bestcon)){
                           bestth=1.0/E; bestg=g; besth=h; bests=s;
                        }
                    }
                }
            }
        
        g=bestg; h=besth; t=g-h+1;
        s=bests;
        calc(&E,&con,1);
        if (1.0/E>maxex) maxex=1.0/E;
        if (con>maxcon) maxcon=con;
    }
    printf("max. ex: %10.6f max. const.: %10.6f\n", maxex,maxcon);
  }
}
