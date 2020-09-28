#include<bits/stdc++.h>
using namespace std;
#define rep(i,n) for(int i=0;i<n;i++)
typedef long long ll;
#define vi vector<ll>
#define all(a) a.begin(),a.end()
typedef pair<int,int> P;
typedef complex<double> C;
typedef long double ld;
const long long mod=1000000007;
const double pi=3.14159265358979;
const double eps=1e-12;
vector<C> a,z,newz,pb;
double comb[20][20];
int n;

C p(vector<C>f, C x){//p(x)を求める
	C ans=0;
	rep(I,f.size())ans+=pow(x,I)*f[I];
	return ans;
}

C pd(vector<C> f,C x){//p'(x)を求める
	C ans=0;
	for(int I=1;I<f.size();I++)ans+=(double)I*pow(x,I-1)*f[I];
	return ans;
}

void disp(C x){//複素数の表示
	cout<<fixed<<setprecision(14)<<x.real()<<" + "<<x.imag()<<" i\n";
}

void input(){//入力　昇冪の順に！
	cin>>n;
	a.resize(n+1);
	rep(i,n+1){
		double re,im;cin>>re>>im;
		a[i]=re+im*1.0i;
	}
	rep(i,n+1)a[i]/=a[n];
	z.resize(n+1);
	newz.resize(n+1);

	comb[0][0]=1;
	rep(i,19){
		rep(j,i+1){
			comb[i+1][j]+=comb[i][j];
			comb[i+1][j+1]+=comb[i][j];
		}
	}

}

C newton(vector<C>f){
	C nx=1000;
	while(abs(p(f,nx))>eps){
		nx-=p(f,nx)/pd(f,nx);
	}
	return nx;
}

int main(){
	input();

	C b=-a[n-1]/(double)n,r=20;
	pb.resize(n+1);
	rep(i,n+1){
		for(int j=0;j<=i;j++){
			pb[j]+=a[i]*comb[i][j]*pow(b,i-j);
		}
	}
	
	rep(i,n)pb[i]=-abs(pb[i]);
	r=newton(pb);
	
	rep(j,n){
		z[j]=b+r*exp(C(1.0i*(pi*2*j/n+3/2/n)));
	}
	int cnt=0,step=0;
	do{
			cnt=0;
			rep(j,n){
			if(abs(p(a,z[j]))<=eps){
					cnt++;
					newz[j]=0;
			}
			else{
					newz[j]=p(a,z[j])/pd(a,z[j]);
					C sum=0;
					rep(k,n){
						if(j==k)continue;
						sum+=1.0/(z[j]-z[k]);
					}
					newz[j]=-newz[j]/(1.0-newz[j]*sum);
					}
			}
			rep(j,n)z[j]+=newz[j];
			step++;
	}while(cnt!=n);
	
	cout<<"step == "<<step<<endl;
	rep(i,n){
			disp(z[i]);
	}

	rep(i,n){
		double err=abs(p(a,z[i]));
		rep(j,n){
			if(i==j)continue;
			err/=abs(z[i]-z[j]);
		}
		cout<<"error : "<<fixed<<setprecision(14)<<err<<endl;
	}
}