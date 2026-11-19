#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <tuple>
#include <numeric>
#include <vector>
using namespace std;
// P(a,b) = i^{popcount(a&b)} X(a) Z(b), bit 0 is x_1.
constexpr int D=128,N=D*D;
int pc(int x){return __builtin_popcount((unsigned)x);}
int amp(int x){int s=0;for(auto t: {7,25,97,42,82,76,52}) if((x&t)==t)s^=1; return s?-1:1;}
int main(){
 array<int,N> w{},aa{},bb{},ii{};
 for(int a=0;a<D;a++){
  array<int,D> z{};for(int x=0;x<D;x++)z[x]=amp(x)*amp(x^a);
  for(int j=1;j<D;j*=2)for(int k=0;k<D;k+=2*j)for(int l=0;l<j;l++){int u=z[k+l],v=z[k+l+j]; z[k+l]=u+v;z[k+l+j]=u-v;}
  for(int b=0;b<D;b++){int p=a*D+b;aa[p]=a;bb[p]=b;ii[p]=pc(a&b);if((ii[p]%2) && z[b]!=0)return 7;w[p]=(ii[p]%2)?0:((ii[p]&2)?-z[b]:z[b]);}
 }
 map<int,int> counts;map<pair<uint64_t,uint64_t>,int> moments,bounds;
 for(int p=1;p<N;p++){
  int a=aa[p],b=bb[p];int dp=D+w[p],dm=D-w[p];if(dp<=0||dm<=0)return 2;
  array<int,257> hp{},hm{};
  for(int m=0;m<N;m++){
   int c=aa[m],d=bb[m]; if((pc(a&d)+pc(b&c))&1){hp[0]++;hm[0]++;continue;}
   int j=p^m;int exp=(ii[m]+ii[p]-ii[j]+2*pc(d&a))&3; if(exp&1)return 3;
   int v=(exp==2?-1:1)*w[j]; int vp=abs(w[m]+v),vm=abs(w[m]-v); if(vp>dp||vm>dm)return 4;hp[vp]++;hm[vm]++;
  }
  for(int e:{1,-1}){
   auto &h=e==1?hp:hm;uint64_t den=e==1?dp:dm;
   if(accumulate(h.begin(),h.end(),0)!=N)return 8;int s=h[den];if(!s||(s&(s-1)))return 5;int eta=7-__builtin_ctz((unsigned)s);counts[eta]++;
   uint64_t purity=0,mom=0; for(uint64_t i=1;i<=den;i++){purity+=h[i]*i*i;mom+=h[i]*i*i*i*i;}
   if(purity!=D*den*den)return 6;
   if(eta==6){uint64_t mden=D*den*den*den*den;auto gcd=std::gcd(mom,mden);moments[{mom/gcd,mden/gcd}]++;
    // Multiply padded target magnitudes by 4 to keep integer numerators.
    int one=2,half=112,quarter=1568;uint64_t sum=0;
    for(int v=den;v>=0;v--){int n=h[v];int use=min(n,one);sum+=4ULL*v*use;n-=use;one-=use;use=min(n,half);sum+=2ULL*v*use;n-=use;half-=use;use=min(n,quarter);sum+=1ULL*v*use;quarter-=use;}
    uint64_t sd=4*D*den;auto sdg=std::gcd(sum,sd);bounds[{sum/sdg,sd/sdg}]++;
   }
  }
 }
 for(auto [e,n]:counts)cout<<"nullity "<<e<<" "<<n<<"\n";
 for(auto [q,n]:moments)cout<<"moment "<<q.first<<"/"<<q.second<<" "<<n<<"\n";
 for(auto [q,n]:bounds)cout<<"bound "<<q.first<<"/"<<q.second<<" "<<n<<"\n";
}
