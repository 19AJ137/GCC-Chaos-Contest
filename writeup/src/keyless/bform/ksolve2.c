/* ksolve2.c — 外側の鍵 k が写像 2 のとき、解いた鍵ストリームの語から写像の定数 Q・A と状態を求める（近似。2026-11-12）
 * 写像 2: a = Q*u; b = a*u; u' = A*u - a*v; v' = C*v + b（C ≈ 1-Q）。出力は v。鍵ストリームの語 = v のバイト (4,2,3,5)。連鎖なし ⇔ バイト 6 > 0xE6。
 * w_n = v_{n+1} - C v_n = Q u_n^2、w_{n+1} = w_n (A - Q v_n)^2 から、連続する 4 語 v_0..v_3 で
 *     s0*sqrt(w1/w0) - s1*sqrt(w2/w1) = Q (v1 - v0)        （w_n = (v_{n+1}-v_n) + Q v_n、s0・s1 = ±1）
 * を Q について解き、A = Q v0 + s0*sqrt(w1/w0)。先の語 v_4・v_5・… を予測して鍵ストリームの語と比べる。
 * v のバイト 6・7（指数と仮数の上位 4 ビット）は全通り試す。
 *   ksolve2 k_words.txt 開始語
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
typedef uint64_t u64; typedef uint32_t u32;
static double u2d(u64 u){double x;memcpy(&x,&u,8);return x;}
static u64 d2u(double x){u64 u;memcpy(&u,&x,8);return u;}
static u32 KW[256]; static int FL[256], HAVE[256];
static double mk(u32 kw,int b6,int b7,int low){ u64 b4=(kw>>24)&0xFF,b2=(kw>>16)&0xFF,b3=(kw>>8)&0xFF,b5=kw&0xFF;
  return u2d((u64)low|(b2<<16)|(b3<<24)|(b4<<32)|(b5<<40)|((u64)b6<<48)|((u64)b7<<56)); }
static u32 kw_of(double x){ u64 u=d2u(x); return (u32)((((u>>32)&0xFF)<<24)|(((u>>16)&0xFF)<<16)|(((u>>24)&0xFF)<<8)|((u>>40)&0xFF)); }
/* 仮数のビット 16〜47 を 32 ビットの数として（上位から b5 b4 b3 b2） */
static u32 mid_of(double x){ u64 u=d2u(x); return (u32)((u>>16)&0xFFFFFFFFu); }
static u32 mid_kw(u32 kw){ u32 b4=(kw>>24)&0xFF,b2=(kw>>16)&0xFF,b3=(kw>>8)&0xFF,b5=kw&0xFF; return (b5<<24)|(b4<<16)|(b3<<8)|b2; }
static int cands(int n,double*out,int*c6,int*c7){ int k=0;
  if(FL[n]==0){ for(int b6=0xE7;b6<=0xFF;b6++){ out[k]=mk(KW[n],b6,0x3F,0x8000); c6[k]=b6; c7[k]=0x3F; k++; } }
  else { for(int b6=0x00;b6<=0x18;b6++){ out[k]=mk(KW[n],b6,0x40,0x8000); c6[k]=b6; c7[k]=0x40; k++; }
         for(int b6=0x60;b6<=0xE6;b6++){ out[k]=mk(KW[n],b6,0x3F,0x8000); c6[k]=b6; c7[k]=0x3F; k++; }
         for(int b6=0x80;b6<=0xD8;b6++){ out[k]=mk(KW[n],b6,0xBF,0x8000); c6[k]=b6; c7[k]=0xBF; k++; } }
  return k; }
static int match(double v,int n,int tol){ if(!(fabs(v)<8.0)) return 0; u32 a=mid_of(v), b=mid_kw(KW[n]); u32 d=a>b?a-b:b-a; if(d>0x80000000u) d=(u32)(0u-d); return d<(u32)tol; }
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); int n,fl; unsigned kw; while(fscanf(f,"%d %d %x",&n,&fl,&kw)==3){ KW[n]=kw; FL[n]=fl; HAVE[n]=1; } fclose(f);
  int n0=atoi(argv[2]); for(int i=0;i<8;i++) if(!HAVE[n0+i]){ fprintf(stderr,"語 %d がない\n",n0+i); return 1; }
  static double C0[512],C1[512],C2[512],C3[512]; static int a6[512],a7[512],b6[512],b7[512],c6[512],c7[512],d6[512],d7[512];
  int k0=cands(n0,C0,a6,a7), k1=cands(n0+1,C1,b6,b7), k2=cands(n0+2,C2,c6,c7), k3=cands(n0+3,C3,d6,d7);
  fprintf(stderr,"候補の数: %d %d %d %d（組 %.2e）\n",k0,k1,k2,k3,(double)k0*k1*k2*k3);
  const int NG=160; const double QLO=0.50, QHI=1.06; long nroot=0, nhit=0;
  for(int i0=0;i0<k0;i0++){ double v0=C0[i0];
   for(int i1=0;i1<k1;i1++){ double v1=C1[i1], d0=v1-v0;
    for(int i2=0;i2<k2;i2++){ double v2=C2[i2], d1=v2-v1;
     for(int i3=0;i3<k3;i3++){ double v3=C3[i3], d2=v3-v2;
      for(int s0=-1;s0<=1;s0+=2) for(int s1=-1;s1<=1;s1+=2){
        double Fp=NAN, Qp=0;
        for(int g=0;g<=NG;g++){ double Q=QLO+(QHI-QLO)*g/NG; double w0=d0+Q*v0, w1=d1+Q*v1, w2=d2+Q*v2; double F;
          if(w0>0&&w1>0&&w2>0) F=s0*sqrt(w1/w0)-s1*sqrt(w2/w1)-Q*d0; else F=NAN;
          if(!isnan(F)&&!isnan(Fp)&&((F<0)!=(Fp<0))){
            double lo=Qp, hi=Q, flo=Fp;
            for(int it=0;it<60;it++){ double m=0.5*(lo+hi); double a0=d0+m*v0,a1=d1+m*v1,a2=d2+m*v2; if(!(a0>0&&a1>0&&a2>0)) break;
              double fm=s0*sqrt(a1/a0)-s1*sqrt(a2/a1)-m*d0; if((fm<0)==(flo<0)){lo=m;flo=fm;} else hi=m; }
            double Qs=0.5*(lo+hi); double w0s=d0+Qs*v0,w1s=d1+Qs*v1,w2s=d2+Qs*v2; if(!(w0s>0&&w1s>0&&w2s>0)) { Fp=F; Qp=Q; continue; }
            double A=Qs*v0+s0*sqrt(w1s/w0s); nroot++;
            if(A>1.2&&A<2.6){
              /* 先の語を予測: w3 = w2 (A - Q v2)^2、v4 = (1-Q) v3 + w3、… */
              double w=w2s, vp=v2, vc=v3; int ok=0; double vv[8];
              for(int j=4;j<8;j++){ w=w*(A-Qs*vp)*(A-Qs*vp); double vn=(1-Qs)*vc+w; vv[j]=vn; if(!match(vn,n0+j,j==4?0x40000:0x4000000)) break; ok++; vp=vc; vc=vn; }
              if(ok>=2){ nhit++; double u0=sqrt(w0s/Qs);
                printf("Q=%.15f A=%.15f |u_%d|=%.12f  v: %.12f %.12f %.12f %.12f  先の語の一致 %d  バイト6/7: %02x%02x %02x%02x %02x%02x %02x%02x  s0=%d s1=%d\n",Qs,A,n0,u0,v0,v1,v2,v3,ok,a6[i0],a7[i0],b6[i1],b7[i1],c6[i2],c7[i2],d6[i3],d7[i3],s0,s1); fflush(stdout); }
            }
          }
          Fp=F; Qp=Q; }
      } } } } }
  fprintf(stderr,"根 %ld 個、一致 %ld 個\n",nroot,nhit); return 0; }
