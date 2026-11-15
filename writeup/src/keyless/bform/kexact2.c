/* kexact2.c — 外側の鍵 k（写像 2）の定数 Q・A・C と状態 (u, v) を、倍精度の値として厳密に求める（2026-11-12）
 * gcc.py の step（mode 2）と同じ順序で計算する:
 *     a = Q*u; b = a*u; nu = A*u - a*v; nv = C*v + b; u = set_low16(nu); v = set_low16(nv)
 * Q・A・C・u・v の下位 16 ビットは鍵で決まる定数 lo（= b0 | b1<<8）。v のバイト 2〜5 は鍵ストリームの語で分かっている。
 * 近似解（ksolve2）のまわりで、lo（65,536 通り）と Q・C・u・A の候補を全部試し、分かっている語が全部合うものを出す。
 *   kexact2 k_words.txt n0 Q A u  b6_0 b7_0 b6_1 b7_1 b6_2 b7_2 b6_3 b7_3  [R]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
typedef uint64_t u64; typedef uint32_t u32;
static inline double u2d(u64 u){double x;memcpy(&x,&u,8);return x;}
static inline u64 d2u(double x){u64 u;memcpy(&u,&x,8);return u;}
static inline double setlow(double x,u64 lo){ return u2d((d2u(x)&~0xFFFFull)|lo); }
static u32 KW[256]; static int FL[256], HAVE[256];
static double mk(u32 kw,int b6,int b7,u64 low){ u64 b4=(kw>>24)&0xFF,b2=(kw>>16)&0xFF,b3=(kw>>8)&0xFF,b5=kw&0xFF;
  return u2d(low|(b2<<16)|(b3<<24)|(b4<<32)|(b5<<40)|((u64)b6<<48)|((u64)b7<<56)); }
static inline u32 kw_of(double x){ u64 u=d2u(x); return (u32)((((u>>32)&0xFF)<<24)|(((u>>16)&0xFF)<<16)|(((u>>24)&0xFF)<<8)|((u>>40)&0xFF)); }
static inline void stepf(double Q,double A,double C,u64 lo,double*u,double*v){
  volatile double a=Q*(*u); volatile double b=a*(*u); volatile double t1=A*(*u); volatile double t2=a*(*v); volatile double nu=t1-t2;
  volatile double t3=C*(*v); volatile double nv=t3+b; *u=setlow(nu,lo); *v=setlow(nv,lo); }
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); int n,fl; unsigned kw; int nmax=0; while(fscanf(f,"%d %d %x",&n,&fl,&kw)==3){ KW[n]=kw; FL[n]=fl; HAVE[n]=1; if(n>nmax)nmax=n; } fclose(f);
  int n0=atoi(argv[2]); double Qa=atof(argv[3]), Aa=atof(argv[4]), Ua=atof(argv[5]); int b6[4],b7[4];
  for(int i=0;i<4;i++){ b6[i]=(int)strtol(argv[6+2*i],0,16); b7[i]=(int)strtol(argv[7+2*i],0,16); }
  long R=argc>14?atol(argv[14]):300; long nsol=0; long st1=0, st2=0;
  for(u64 lo=0;lo<65536;lo++){
    double V0=mk(KW[n0],b6[0],b7[0],lo), V1=mk(KW[n0+1],b6[1],b7[1],lo), V2=mk(KW[n0+2],b6[2],b7[2],lo), V3=mk(KW[n0+3],b6[3],b7[3],lo);
    u64 q0=(d2u(Qa)&~0xFFFFull)|lo;
    for(long dq=-R;dq<=R;dq++){ double Q=u2d(q0+(u64)(dq*0x10000));
      double c0=1.0-Q; u64 cu=(d2u(c0)&~0xFFFFull)|lo;
      for(long dc=-6;dc<=6;dc++){ double C=u2d(cu+(u64)(dc*0x10000));
        /* u の候補: V1 = setlow(C*V0 + Q*u*u) から */
        double w=(V1-C*V0)/Q; if(!(w>0)) continue; double ue=sqrt(w); if(fabs(ue-Ua)>1e-6) continue;
        u64 uu=(d2u(ue)&~0xFFFFull)|lo;
        for(long du=-3;du<=3;du++){ double U=u2d(uu+(u64)(du*0x10000));
          volatile double a=Q*U; volatile double b=a*U; volatile double t3=C*V0; volatile double nv=t3+b;
          if(setlow(nv,lo)!=V1) continue; st1++;
          u64 au=(d2u(Aa)&~0xFFFFull)|lo;
          for(long da=-R;da<=R;da++){ double A=u2d(au+(u64)(da*0x10000));
            double u=U, v=V0; stepf(Q,A,C,lo,&u,&v);        /* → (U1, V1) */
            stepf(Q,A,C,lo,&u,&v); if(v!=V2) continue;      /* → V2 */
            stepf(Q,A,C,lo,&u,&v); if(v!=V3) continue; st2++;
            /* 分かっている語を全部照合 */
            int ok=0,bad=0; for(int m=n0+4;m<=nmax;m++){ stepf(Q,A,C,lo,&u,&v); if(HAVE[m]){ if(kw_of(v)==KW[m]) ok++; else bad++; } }
            if(bad<=3){ nsol++; printf("lo=%04llx Q=%.17g (%016llx) A=%.17g (%016llx) C=%.17g (%016llx) u_%d=%.17g (%016llx) v_%d=%.17g (%016llx)  一致 %d 不一致 %d\n",(unsigned long long)lo,Q,(unsigned long long)d2u(Q),A,(unsigned long long)d2u(A),C,(unsigned long long)d2u(C),n0,U,(unsigned long long)d2u(U),n0,V0,(unsigned long long)d2u(V0),ok,bad); fflush(stdout); }
          } } } } }
  fprintf(stderr,"1 段目の通過 %ld、3 語の通過 %ld、解 %ld\n",st1,st2,nsol); return 0; }
