// mainprogram1.c and mainprogram2.c combined
// implemented using ARB Ball arithmetic
#include "stdbool.h"
#include "arb.h"

long my_floor(arb_t x, long prec)
{
  static arf_t lb,ub;
  static bool init=false;
  if(!init)
    {
      init=true;
      arf_init(lb);arf_init(ub);
    }
  arb_get_lbound_arf(lb,x,prec);
  arb_get_ubound_arf(ub,x,prec);
  long llb=arf_get_si(lb,ARF_RND_DOWN);
  long lub=arf_get_si(ub,ARF_RND_DOWN);
  if(llb!=lub)
    {
      printf("Fatal error in my_floor. Exiting.\n");
      exit(0);
    }
  return llb;
}

void newdel(arb_t res, long k, long r, arb_t del, long prec)
{
  long tkr,jj,j;
  static arb_t y,p,tmp1,tmp2,tmp3,tmp4;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(y);arb_init(p);arb_init(tmp1);arb_init(tmp2);
      arb_init(tmp3);arb_init(tmp4);
    }
  if((r<4) && (r> k))
    {
      arb_mul_2exp_si(res,del,1);
      return;
    }
  tkr=2*k*r;
  arb_mul_2exp_si(tmp1,del,1);
  arb_sub_ui(y,tmp1,(k-r)*(k-r+1),prec); // y=2del-(k-r)(k-r+1)
  //printf("y=");arb_printd(y,20);printf("\n");
  if(arb_is_negative(y))
    {
      arb_set_ui(tmp1,k+1);
      arb_inv(tmp2,tmp1,prec); // 1/(k+1)
      arb_add_ui(tmp1,y,tkr,prec);
      arb_set_ui(tmp3,k*2);
      arb_div(tmp4,tmp3,tmp1,prec); // 2k/(tkr+y)
      arb_sub(tmp3,tmp4,tmp2,prec);
      if(!arb_is_positive(tmp3)) // 2k/(tkr+y) <= 1/(k+1)
	{
	  arb_mul_2exp_si(res,del,1);
	  return;
	}
    }
  arb_set_ui(tmp1,9*r);
  arb_div_ui(tmp2,tmp1,10,prec);
  arb_mul_2exp_si(tmp1,y,2); // 4y
  arb_add_ui(tmp3,tmp1,1,prec); // 4y+1
  arb_sqrt(tmp1,tmp3,prec); // sqrt(.)
  arb_add_ui(tmp3,tmp1,3,prec); // 3+sqrt(.)
  arb_mul_2exp_si(tmp3,tmp3,-1); // 0.5(3+sqrt(.))
  arb_min(tmp1,tmp3,tmp2,prec);
  j=my_floor(tmp1,prec);
  //printf("j=%ld\n",j);
  arb_set_ui(tmp1,r);
  arb_inv(p,tmp1,prec); // 1/r
  arb_mul_2exp_si(tmp4,p,-1); // 0.5/r
  for(jj=j-1;jj>=1;jj--)
    {
      arb_div_ui(tmp1,y,tkr-2*r*jj,prec);
      arb_sub_ui(tmp2,tmp1,1,prec);
      arb_neg(tmp2,tmp2);
      arb_mul(tmp1,tmp2,p,prec);
      arb_mul_2exp_si(tmp1,tmp1,-1);
      arb_add(p,tmp4,tmp1,prec);
    }
  //printf("p=");arb_printd(p,20);printf("\n");
  arb_sub_ui(tmp1,y,tkr,prec); // (y-tkr)
  arb_mul(tmp2,tmp1,p,prec); // (y-tkr)p
  arb_mul_2exp_si(tmp2,tmp2,-1); // (y-tkr)p/2
  arb_sub(tmp1,del,tmp2,prec); // return del+(.)
  arb_sub_ui(res,tmp1,k,prec);
  return;
}



int main(int argc, char **argv)
{
  printf("Command:- ");
  for(int i=0;i<argc;i++)
    printf("%s ", argv[i]);
  printf("\n");
  if(argc!=5)
    {
      printf("Usage:- %s <k0> <k1> <1 or 2> <prec>\n",argv[0]);
      return 0;
    }
  bool prog1;
  if(atol(argv[3])==1)
    prog1=true;
  else
    prog1=false;
  long prec=atol(argv[4]);
  long j,k,k0,k1,r,r0,r1,n,bestr,s;
  arb_t kk,logk,del0,del1,bestdel,goal,maxs,eta,om;
  arb_t logH,logW,logC,k3,theta,thetamax;
  arb_t tmp1,tmp2,tmp3,tmp4;
  arb_init(kk);arb_init(logk);arb_init(del0);arb_init(del1);arb_init(bestdel);
  arb_init(goal);arb_init(maxs);arb_init(eta);arb_init(om);arb_init(logH);
  arb_init(logW);arb_init(logC);arb_init(k3);arb_init(theta);arb_init(thetamax);
  arb_init(tmp1);arb_init(tmp2);arb_init(tmp3);arb_init(tmp4);
  k0=atol(argv[1]);
  k1=atol(argv[2]);
  for(k=k0;k<=k1;k++)
    {
      arb_log_ui(logk,k,prec);
      arb_set_ui(tmp1,k*k*k);arb_mul(k3,tmp1,logk,prec);
      arb_set_d(om,0.5);
      arb_set_d(tmp3,1.5);
      for(j=1;j<=10;j++)
	{
	  arb_mul_ui(tmp1,k3,18,prec); // 18k3
	  arb_div(tmp2,tmp1,om,prec); // 18*k3/om
	  arb_log(tmp1,tmp2,prec); // log(18*k3/om)
	  arb_sub(tmp2,tmp1,tmp3,prec); // log(18*k3/om)-1.5
	  arb_div(om,tmp3,tmp2,prec); // 1.5/(log(.)-1.5)
	}
      arb_add_ui(eta,om,1,prec);
      arb_div(tmp1,tmp3,om,prec);
      arb_add(tmp2,tmp1,tmp3,prec); // 1.5+1.5/om
      arb_mul_ui(tmp1,k3,18,prec);
      arb_div(tmp3,tmp1,om,prec);
      arb_log(tmp1,tmp3,prec);
      arb_max(tmp3,tmp1,tmp2,prec);
      arb_mul_ui(logW,tmp3,k+1,prec);
      if(prog1) // 0.5 k^2 (1-1/k)
	{
	  arb_set_ui(tmp1,k);
	  arb_inv(tmp2,tmp1,prec); // 1/k
	  arb_sub_ui(tmp1,tmp2,1,prec); // 1/k-1
	  arb_neg(tmp1,tmp1); // 1-1/k
	  arb_mul_ui(del0,tmp1,k*k,prec); // k^2(1-1/k)
	  arb_mul_2exp_si(del0,del0,-1); // 0.5 k^2 (1-1/k)
	}
      else
	{
	  arb_set_ui(tmp1,4*k*k);
	  arb_div_ui(del0,tmp1,10,prec); // 0.4 k^2
	}
      arb_set_ui(tmp1,k*k);
      arb_div_ui(goal,tmp1,1000,prec);
      arb_log(tmp1,eta,prec);
      arb_mul_ui(tmp2,tmp1,k*k-4*k,prec);
      arb_mul_ui(tmp1,logk,3*k,prec);
      arb_add(logH,tmp1,tmp2,prec);
      arb_mul_ui(logC,logk,k,prec);
      long n0=1;
      if(!prog1)
	n0=ceil(0.1247*k);
      for(n=n0;;n++)
	{
	  //if((n%1000)==0){printf("\nn=%ld del0=",n);arb_printd(del0,20);printf("\n");fflush(stdout);}
	  arb_mul_2exp_si(tmp1,del0,1); // 2del0
	  arb_sub_ui(tmp2,tmp1,k*k+k,prec); // 2del0-k*k+k
	  arb_neg(tmp2,tmp2);  // k*k-k+2del0
	  arb_sqrt(tmp1,tmp2,prec); // sqrt(.)
	  arb_set_d(tmp2,0.5);
	  arb_add(tmp3,tmp1,tmp2,prec); // 0.5+sqrt(.)
	  r0=my_floor(tmp3,prec);
	  r0-=2;
	  r1=r0+4;
	  arb_set_ui(bestdel,k*k);
	  bestr=-1;
	  for(r=r0;r<=r1;r++)
	    {
	      newdel(del1,k,r,del0,prec);
	      //printf("newdel(%ld,%ld,",k,r);arb_printd(del0,10);printf(") returned ");arb_printd(del1,10);printf("\n"); 
	      arb_sub(tmp1,del1,bestdel,prec);
	      if(arb_is_negative(tmp1))
		{
		  arb_set(bestdel,del1);
		  bestr=r;
		}
	    }
	  arb_set(del1,bestdel);
	  r=bestr;
	  arb_sub(tmp1,del0,del1,prec);
	  if(arb_is_negative(tmp1)&&(r<r0)) exit(-1);
	  arb_mul(tmp2,tmp1,logW,prec);
	  arb_log(tmp3,eta,prec);
	  arb_mul_ui(tmp4,tmp3,4*k*n,prec);
	  arb_add(tmp3,tmp4,logH,prec);
	  arb_max(tmp4,tmp3,tmp2,prec);
	  arb_add(logC,logC,tmp4,prec);
	  arb_sub(tmp2,del1,goal,prec);
	  if(arb_is_negative(tmp2))
	    {
	      printf("%ld %ld: s=",k,n);
	      arb_sub(tmp2,del0,goal,prec);
	      arb_div(tmp3,tmp2,tmp1,prec);
	      arb_add_ui(tmp1,tmp3,n,prec);
	      arb_mul_ui(tmp2,tmp1,k,prec);
	      arb_add_ui(tmp1,tmp2,1,prec);
	      s=my_floor(tmp1,prec);
	      arb_div(theta,logC,k3,prec);
	      arb_set_ui(tmp1,s);
	      arb_div_ui(tmp2,tmp1,k*k,prec);
	      arb_printd(tmp2,20);printf(" eta=");arb_printd(eta,20);
	      printf(" theta=");arb_printd(theta,20);printf("\n");
	      arb_max(maxs,maxs,tmp2,prec);
	      arb_max(thetamax,thetamax,theta,prec);
	      break;
	    }
	  arb_set(del0,del1);
	}
    }
  printf("maxs = ");arb_printd(maxs,20);printf("\n");
  printf("maxtheta = ");arb_printd(thetamax,20);printf("\n");
  return 0;
}
