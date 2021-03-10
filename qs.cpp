#include <bits/stdc++.h>

// TODO: debug this (works for several cases but seems to fail for even semiprimes)

using namespace std;

typedef long long ll; 
// consider using __int_128 or some bigint class (https://gist.github.com/vinayak-garg/4007974) for future iterations.
// Could use GMP (https://gmplib.org/) but not supported by online competitive programming judges. Also numpy (for a python version)
// could be good but is not supported by DMOJ

// Sources used:
// https://www.cc.gatech.edu/~rpeng/CS6550_S21/Feb15QS.pdf, https://www.cs.virginia.edu/crab/QFS_Simple.pdf, https://en.wikipedia.org/wiki/Quadratic_sieve


ll n, m; // assume <= 2e18 for now

// parameters to tweak
const int prime_bound = 40; // smoothness bound (i.e. largest prime in the factorization of each number)
const int max_prime_freq = 40; // upper bound for number of primes used
const int sieve_bound = 40; // look at factorizations of x^2 - n where x \in [m + 1, m + sieve_bound] where m = floor(sqrt(n))
const int max_basis_size = 40; // upper bound for the number of smooth elements
const int max_exp = 30; // max exponent of a prime number

vector<int> primes;
vector<vector<ll>> sqrts;
int num_primes = 0;
int is_prime[prime_bound];
vector<vector<ll>> all_powers;

void process_primes() {
    for (int i = 2; i < prime_bound; i++) {
        is_prime[i] = 1;
    }
    for (int i = 2; i < prime_bound; i++) {
        if (is_prime[i]) {
            primes.push_back(i);
            num_primes++;
            for (int j = 2 * i; j < prime_bound; j += i) {
                is_prime[j] = 0;
            }
        }
    }
    for (int p: primes) {
        vector<ll> powers = {1};
        for (int i = 1; i < max_exp; i++) {
            powers.push_back((powers[i - 1] * p) % n);
        }
        all_powers.push_back(powers);
    }
}

void find_sqrts() {
    sqrts.resize(num_primes);
    for (int i = 0; i < num_primes; i++) {
        int p = primes[i];
        int nn = n % p;
        // find sqrts of n mod p
        for (int j = 0; j < p; j++) {
            if ((j * j) % p == nn) {
                sqrts[i].push_back(j);
            }
        }
    }
}

int basis_size = 0;

vector<ll> sieved_sqrts; // sqrts of x^2
vector<vector<int>> sieved_alts; // factors of x^2 - n
vector<bitset<max_basis_size>> basis;
vector<ll> prime_divs;

void sieve() {
    vector<ll> elts; // elements of the form (x + m)^2 - n (with x \in [1, m + sieve_bound]) is more convenient to implement
    vector<vector<int>> prime_powers(sieve_bound, vector<int>(num_primes));
    int l = m + 1;
    int r = m + sieve_bound;
    for (ll i = l; i <= r; i++) {
        elts.push_back((i * i - n));
    }
    vector<ll> elts0 = elts; // keep track of original values (for debugging)
    for (int i = 0; i < num_primes; i++) {
        int p = primes[i];
        for (int sq: sqrts[i]) {
            // want p | x^2 - n
            // i.e. x^2 \equiv n mod p
            // solutions of x are (s0 + pk), (s1 + pk), where s0 and s1 are square roots of n mod p
            int min_idx = l - (l % p) + sq; // min idx = sq + pk
            if (min_idx < l) min_idx += p;
            assert(min_idx >= l);
            for (int j = min_idx; j <= r; j += p) {
                int idx = j - l;
                assert(idx >= 0 and idx <= r - l);
                while (elts[idx] % p == 0) {
                    elts[idx] /= p;
                    prime_powers[idx][i]++;
                    if (elts[idx] == 1) {
                        // j^2 - n is smooth, add to basis
                        prime_divs.push_back(p);
                    }
                }
            }
        }
    }

    basis.resize(num_primes);
    int cur_basis = 0;
    for (int i = l; i <= r; i++) {
        int idx = i - l;
        if (elts[idx] == 1) {
            // i^2 is smooth, add to basis
            sieved_sqrts.push_back(i);
            sieved_alts.emplace_back();
            sieved_alts.back().resize(num_primes);
            for (int j = 0; j < num_primes; j++) {
                int pw = prime_powers[idx][j];
                sieved_alts.back()[j] = pw;
            }
            for (int j = 0; j < num_primes; j++) {
                basis[j][cur_basis] = prime_powers[idx][j] & 1;
            }
            cur_basis++;
        }
    }
    basis_size = (int)sieved_sqrts.size();
}
vector<bitset<max_basis_size>> null_space;

void print_basis() {
    for (int row = 0; row < num_primes; row++) {
        cout << basis[row] << "\n";
    }
    cout << "\n";
}

void find_null_space() {
    // Find solutions to Ax = 0
    // Size of A is num prime factors x basis size
    vector<int> where(basis_size, -1);
    // print_basis();
    for (int row = 0, col = 0; col < basis_size and row < num_primes; col++) {
        // cout << "col: " << col << ", row: " << row << "\n";
        for (int i = row; i < num_primes; i++) {
            if (basis[i][col]) {
                swap(basis[i], basis[row]);
                where[col] = row;
                break;
            }
        }
        if (basis[row][col]) {
            for (int i = 0; i < row; i++) {
                if (basis[i][col]) basis[i] ^= basis[row];
            }
            for (int i = row + 1; i < num_primes; i++) {
                if (basis[i][col]) basis[i] ^= basis[row];
            }
            row++;
        }
        // cout << "where: " << where[col] << "\n";
        // print_basis();
    }
    // print_basis();
    // cout << "creating null space\n";
    for (int col = 0; col < basis_size; col++) {
        if (where[col] == -1) {
            bitset<max_basis_size> bs;
            for (int c2 = 0; c2 < col; c2++) {
                if (~where[c2]) {
                    int row = where[c2];
                    bs[c2] = basis[row][col];
                } 
            }
            bs[col] = 1;
            null_space.push_back(bs);
            // cout << "col: " << col << ", bs: " << bs << "\n";
        }
    }
}

ll p = -1;
ll q = -1;

void try_coeffs() {
    for (int it = 0; it < (int)null_space.size(); it++) {
        bitset<max_basis_size>& bs = null_space[it];
        // cout << "bs: " << bs << "\n";
        ll a = 1;
        ll b = 1;
        vector<int> b_vec(num_primes);
        for (int i = 0; i < basis_size; i++) {
            if (bs[i]) {
                a *= sieved_sqrts[i];
                for (int j = 0; j < num_primes; j++) {
                    int pw = sieved_alts[i][j];
                    b_vec[j] += pw;
                }
                a %= n;
            }
        }
        for (int i = 0; i < num_primes; i++) {
            b_vec[i] /= 2;
            int pw = b_vec[i];
            b *= all_powers[i][pw];
            b %= n;
        }
        // a^2 = b^2 mod n
        ll pp = (a + b) % n;
        ll qq = (a - b) % n;
        if (qq < 0) qq += n;
        
        if (abs(a) == abs(b)) continue;
        // cout << "a: " << a << ", b: " << b << "\n";
        if (pp > 1 and pp < n and qq > 1 and qq < n and pp != qq and pp * qq == n) {
            // try all four preimages
            ll pp2 = n - pp;
            ll qq2 = n - qq;
            pair<ll, ll> opt = {-1, -1};
            for (int a1: (vector<ll>){pp, pp2}) {
                for (int a2: (vector<ll>){qq, qq2}) {
                    if (a1 * a2 == n) {
                        opt = {a1, a2};
                    }
                }
            }
            if (~opt.first) {
                p = min(opt.first, opt.second), q = min(opt.first, opt.second);
            }
            p = min(pp, qq), q = max(pp, qq);
            return;
        }
    }
}

int main() {
    // freopen("A.in", "r", stdin);
    // freopen("A.out", "w", stdout);
    cin >> n;
    m = (ll)sqrt(n);

    // cout << "m: " << m << "\n";

    process_primes();

    find_sqrts(); // find sqrts of n mod p for each p

    sieve(); // find prime factorizations of f(x) for all x in [m + 1, m + sieve_bound]

    find_null_space(); // find solutions to Ax = 0

    try_coeffs(); // iterate over a, b such that a^2 = b^2 mod n

    if (~p) {
        cout << n << " " << p << " " << q << "\n";
    } else {
        cout << "No factorization found!\n";
    }
    
}