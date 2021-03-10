#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

int t;

const int maxp = 89;

int isprime[maxp];

int main() {
    ios::sync_with_stdio(0);
    cin.tie(0);
    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

    for (int i = 2; i < maxp; i++) {
        isprime[i] = 1;
    }
    vector<ll> primes;
    for (int i  = 2; i < maxp; i++) {
        if (isprime[i]) {
            for (int j = 2 * i; j < maxp; j += i) {
                isprime[j] = 0;
            }
            if (i != 2) primes.push_back(i);
        }
    }
    int num_primes = primes.size();
    int i1 = rng() % num_primes;
    int i2;
    while (true) {
        i2 = rng() % num_primes;
        if (i2 != i1) break;
    }
    ll p = primes[i1];
    ll q = primes[i2];
    ll n = p * q;
    cout << n << "\n";
} 