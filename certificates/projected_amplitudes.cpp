#include <array>
#include <iostream>
#include <map>
#include <numeric>
#include <cassert>
#include <cstdint>
using namespace std;
using I=long long; using Fr=pair<I,I>;
Fr reduce(I p,I q){I g=gcd(p,q);return {p/g,q/g};}
int parity(unsigned a){return __builtin_parity(a);} int pop(unsigned a){return __builtin_popcount(a);}
int main(){constexpr int D=128;array<int,D> s{}; int terms[]={7,25,97,42,82,76,52};
 for(int x=0;x<D;x++){int count=0;for(int t:terms) count^=( (t&x)==t );s[x]=1-2*count;}
 map<int,int> nulls;map<Fr,int> moments,caps;Fr best{0,1};int branch=0;
 for(int a=0;a<D;a++)for(int b=0;b<D;b++)if(a||b)for(int e:{-1,1}){
  array<int,D> re{},im{};int norm=0,phase=pop(a&b)%4;
  for(int x=0;x<D;x++){int k=e*s[x^a]*(1-2*parity(b&(x^a)));re[x]=s[x]+(phase==0?k:phase==2?-k:0);im[x]=(phase==1?k:phase ==3?-k:0);norm+=re[x]*re[x]+im[x]*im[x];}
  assert(norm>0);array<int,1025> hist{};I fourth=0,squared=0;int stab=0;
  for(int d=0;d<D;d++){
   array<int,D> R{},J{};for(int x=0;x<D;x++){int y=x^d;R[x]=re[y]*re[x]+im[y]*im[x];J[x]=re[y]*im[x]-im[y]*re[x];}
   for(int len=1;len<D;len*=2)for(int x=0;x<D;x+=2*len)for(int z=x;z<x+len;z++){int q=R[z],r=R[z+len],u=J[z],v=J[z+len];R[z ]=q+r;R[z+len]=q-r;J[z]=u+v;J[z+len]=u-v;}
   for(int h=0;h<D;h++){bool odd=parity(d&h);assert((odd?R[h]:J[h])==0);I val=abs(odd?J[h]:R[h]);assert(val<=norm);hist[val ]++;stab+=(val==norm);squared+=val*val;fourth+=val*val*val*val;}
  }
  assert(accumulate(hist.begin(),hist.end(),0)==D*D);assert(squared==I(D)*norm*norm);assert(stab&&!(stab&(stab-1)));int n=7-__builtin_ctz((unsigned)stab);nulls[n]++;branch++; cout<<"BR "<<a<<" "<<b<<" "<<e<<" "<<n<<" "<<norm/2<<"\n";
  if(n==6){moments[reduce(fourth,I(D)*norm*norm*norm*norm)]++;
    int cum=0;I score=0;for(int val=norm;val>0;--val){int next=cum+hist[val];auto mass=[&](int lo,int hi){return max(0,min( next,hi)-max(cum,lo));};score+=I(val)*(4*mass(0,2)+2*mass(2,114)+mass(114,1682));cum=next;}
    Fr cap=reduce(score,4LL*D*norm);caps[cap]++; auto mm=reduce(fourth,I(D)*norm*norm*norm*norm); cout<<"CERT "<<a<<" "<<b<<" "<<e<<" "<<mm.first<<"/"<<mm.second<<" "<<cap.first<<"/"<<cap.second<<"\n";if(cap.first*best.second>best.first*cap.second)best=cap;
  }
 }
 cout<<"branches "<<branch<<"\n";for(auto [n,c]:nulls)cout<<"nullity "<<n<<" "<<c<<"\n";for(auto [f,c]:moments)cout<<"P2 " <<f.first<<"/"<<f.second<<" "<<c<<"\n";for(auto [f,c]:caps)cout<<"spectrum bound "<<f.first<<"/"<<f.second<<" "<<c<<"\n";cout<<"max spectrum bound "<<best.first<<"/"<<best.second<<"\n";
 assert(branch==32766 && best==Fr(121,128));
}
