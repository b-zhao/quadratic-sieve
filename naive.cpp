#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

// can replace this later with pollards-rho

int main() {
    ios::sync_with_stdio(0);
    cin.tie(0);
    ll n;
    cin >> n;
    for (ll i = 2; i * i < n; i++) {
        if (n % i == 0) {
            cout << n << " " << i << " " << (n/i) << "\n";
            return 0;
        }
    }
} 