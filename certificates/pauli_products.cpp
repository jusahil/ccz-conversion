#include <iostream>
#include <vector>
#include <map>
#include <numeric>
#include <array>
#include <cstdint>
#include <cassert>
using namespace std;
using LL=long long; using U=uint64_t; using Fr=pair<LL,LL>;
Fr frac(LL a,LL b){LL g=gcd(a,b);return {a/g,b/g};}
int pc(int x){return __builtin_popcount((unsigned)x);} int par(int x){return pc(x)&1;}
U fourth(int x){U y=abs(x);return y*y*y*y;}
int main(){constexpr int D=128,N=D*D; vector<int>w(N),phase(N),a(N),b(N),s(D);
 int terms[]={7,25,97,42,82,76,52};
 for(int x=0;x<D;x++){int p=0;for(int t:terms)p^=((x&t)==t);s[x]=1-2*p;}
 for(int i=0;i<N;i++){a[i]=i/D;b[i]=i%D;phase[i]=pc(a[i]&b[i]);int z=0;for(int x=0;x<D;x++)z+=s[x]*s[x^a[i]]*(1-2*par(b[i]& x));if(phase[i]%2)assert(z==0);w[i]=phase[i]%2?0:((phase[i]%4==2)?-z:z);}
 map<int,int> cnt;map<Fr,int> moments, caps; Fr best={0,1};int bestP=0,bestS=0;int all=0;
 for(int p=1;p<N;p++)for(int e:{-1,1}){int den=D+e*w[p];assert(den>0);array<int,257> hist{};U m4=0;int stab=0;
 for(int m=0;m<N;m++){if(par((a[m]&b[p])^(a[p]&b[m]))){hist[0]++;continue;}int j=m^p;int k=(phase[m]+phase[p]-phase[j]+2*pc(b[m]&a[p])) &3;assert(k%2==0);int v=abs(w[m]+e*(k==0?1:-1)*w[j]);assert(v<=den);hist[v]++;m4+=fourth(v);stab+=(v==den);}
 assert(accumulate(hist.begin(),hist.end(),0)==N); U purity=0; for(int v=0;v<=den;v++) purity+=U(hist[v])*v*v; assert(purity==U(D)*den*den); assert(stab>0 && (stab&(stab-1))==0);int nul=7-__builtin_ctz((unsigned)stab);cnt[nul]++;all++; cout<<"BR "<<a[p]<<" "<<b[p]<<" "<<e<<" "<<nul<<" "<<den<<"\n";
 if(nul==6){Fr mo=frac(m4,U(D)*fourth(den));moments[mo]++;assert(mo!=frac(121,1024));
 int remaining[3]={2,112,1568},weights[3]={4,2,1},g=0;LL score=0;
 for(int v=den;v>=0 && g<3;v--){int c=hist[v];while(c && g<3){int take=min(c,remaining[g]);score+=LL(take)*v*weights[g]; remaining[g]-=take;c-=take;if(!remaining[g])g++;}}
 Fr bound=frac(score,4LL*D*den);caps[bound]++; cout<<"CERT "<<a[p]<<" "<<b[p]<<" "<<e<<" "<<mo.first<<"/"<<mo.second<<" "<<bound.first<<"/"<<bound.second<<"\n";if(bound.first*best.second>best.first*bound.second){best=bound;bestP=p;bestS= e;}
 }}
 cout<<"branches "<<all<<"\n";for(auto [k,v]:cnt)cout<<"nullity "<<k<<" "<<v<<"\n";for(auto [f,v]:moments)cout<<"P2 "<<f. first<<"/"<<f.second<<" "<<v<<"\n";for(auto [f,v]:caps)cout<<"spectrum bound "<<f.first<<"/"<<f.second<<" "<<v<<"\n"; cout<<"max spectrum bound "<<best.first<<"/"<<best.second<<" witness (a,b,sign) "<<bestP/D<<","<<bestP%D<<","<<bestS <<"\n";
 assert(all==32766 && cnt[6]==31644 && best.first<best.second);
}
