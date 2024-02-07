#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define max(x,y) (((x)>(y))?(x):(y))
#define min(x,y) (((x)>(y))?(y):(x))

double newdel(double k, double r, double del)
{
    double y,p,tkr;
    long j,jj;
    if ((r<4.0) && (r>k)) return (2.0*del);
    tkr=2.0*k*r; y=2.0*del-(k-r)*(k-r+1.0);
    if ((y<0.0) && (2.0*k/(tkr+y))<= 1.0/(k+1.0))
       return (del*2.0);

    j=floor(min(0.5*(3.0+sqrt(4.0*y+1.0)), 9.0*r/10.0));
    p=1.0/r;
    for (jj=j-1; jj>=1; jj--) {
        p=0.5/r+0.5*(1.0-y/(tkr-2.0*r*jj))*p;
    }
    return(del-k+0.5*p*(tkr-y));
}

int main()
{
    long j,k,k0,k1,r,r0,r1,n,bestr,s;
    double kk,logk,del0,del1,bestdel,goal,maxs,eta,om;
    double logH,logW,logC,k3,theta,thetamax;
    printf("enter k range :"); scanf("%ld %ld", &k0,&k1);
    maxs=0.0; thetamax=0.0;
    for (k=k0; k<=k1; k++) {
        kk=(double)k;
        logk=log(kk); k3=kk*kk*kk*logk;
        om=0.5;
        for (j=1;j<=10;j++) om=1.5/(log(18.0*k3/om)-1.5);
        eta=1.0+om;
        logW=(kk+1.0)*max(1.5+1.5/om,log(18.0/om*k3));
        del0=0.5*kk*kk*(1.0-1.0/kk);
        goal=0.001*kk*kk;
        logH=3.0*kk*logk+(kk*kk-4.0*kk)*log(eta);
        logC=kk*logk;
        for (n=1; ;n++) {
            r0=(long)(sqrt(kk*kk+kk-2.0*del0)+0.5)-2;
            r1=r0+4;
            bestdel=kk*kk; bestr=-1;
            for (r=r0;r<=r1;r++) {
                del1=newdel(kk,(double)r,del0);
                if (del1<bestdel) { bestdel=del1; bestr=r;}
            }
            del1=bestdel; r=bestr;
            if ((del1>=del0) && (r<r0)) exit(-1);
            logC +=max(logH+4.0*kk*n*log(eta),logW*(del0-del1));
            if (del1<=goal) {
                s=(long)((n+(del0-goal)/(del0-del1))*kk+1);
                theta=logC/k3;
                printf("%4d: s=%8.6f k^2 eta=%9.7f theta=%10.8f\n",
                       k,s/kk/kk,eta,theta);
                if ((s/kk/kk)>maxs) maxs=s/kk/kk;
                if (theta>thetamax) thetamax=theta;
                break;
            }
            del0=del1;
        }
    }
    printf("\n max s=%9.6fk^2 maxtheta=%10.8f\n",maxs,thetamax);
    system("pause");
}

