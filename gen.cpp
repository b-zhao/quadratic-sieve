#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

int t;

// const int maxp = 5e7;

// int isprime[maxp];

// source: https://math.stackexchange.com/questions/424238/random-primes-between-4000000000-and-4294967291-c

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
int isprime(unsigned long long n) {
  /*if((n&1)==0) return n==2;*/
  if(n%3==0) return n==3;
  /*if(n<25) return n>1;*/
  unsigned long long p = 5;
  while (p*p <= n) {
    if (n%p==0) return 0;
    p += 2;
    if (n%p==0) return 0;
    p += 4;
  }
  return 1;
}

unsigned long long rand_prime(int lower, int upper) {
  unsigned long long spread = upper - lower + 1;
  while(1) {
    unsigned long long p = 1 | (rng() % spread + lower);
    if (isprime(p)) return p;
  }
}

int main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    // for (int i = 2; i < maxp; i++) {
    //     isprime[i] = 1;
    // }
    // vector<ll> primes;
    // for (int i  = 2; i < maxp; i++) {
    //     if (isprime[i]) {
    //         for (int j = 2 * i; j < maxp; j += i) {
    //             isprime[j] = 0;
    //         }
    //         if (i != 2) primes.push_back(i);
    //     }
    // }
    // int num_primes = primes.size();
    // int i1 = rng() % num_primes;
    // int i2;
    // while (true) {
    //     i2 = rng() % num_primes;
    //     if (i2 != i1) break;
    // }
    int lower = 5e7;
    int upper = 1e9;
    // int lower =1e4;
    // int upper = 3e5;
    ll p = rand_prime(lower, upper);
    ll q = rand_prime(lower, upper);
    ll n = p * q;
    cout << n << "\n";
} 