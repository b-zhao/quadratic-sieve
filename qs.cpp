#include <bits/stdc++.h>


// To focus on: How to efficiently get solutions from the basis vectors of the solutions to Ax = 0? (right now I am just
// iterating through the basis vectors of the null space, while I could potentially be using information from the entire
// null space). Maybe some kind of random sampling of subsets would be good?


using namespace std;


typedef int int2; 
#define int __int128 // this prevents overflow for p, q larger than 1e6
typedef __int128 ll; 

// consider using __int_128 or some bigint class (https://gist.github.com/vinayak-garg/4007974) for future iterations.
// Could use GMP (https://gmplib.org/) but not supported by online competitive programming judges. Also numpy (for a python version)
// could be good but is not supported by DMOJ

// Sources used:
// https://www.cc.gatech.edu/~rpeng/CS6550_S21/Feb15QS.pdf, https://www.cs.virginia.edu/crab/QFS_Simple.pdf, https://en.wikipedia.org/wiki/Quadratic_sieve

// 599607257847593 = 15418259 * 38889427

ll n, m; // assume <= 2e18 for now

// parameters to tweak
const int prime_bound = 3000; // smoothness bound (i.e. largest prime in the factorization of each number)
const int sieve_bound = 50000; // look at factorizations of x^2 - n where x \in [m + 1, m + sieve_bound] where m = floor(sqrt(n))
const int max_basis_size = 50000; // upper bound for the number of smooth elements
const int max_exp = 400; // max exponent of a prime number

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
vector<vector<int>> sieved_alts; // prime factorizations of x^2 - n
vector<bitset<max_basis_size>> basis;

void sieve() {
    vector<ll> elts; // elements of the form x^2 - n (with x \in [m + 1, m + sieve_bound])
    vector<vector<int>> prime_powers(sieve_bound, vector<int>(num_primes));
    int l = m + 1;
    int r = m + sieve_bound;
    for (ll i = l; i <= r; i++) {
        elts.push_back(i * i - n);
    }
    vector<ll> elts0 = elts; // keep track of original values (for debugging)
    for (int i = 0; i < num_primes; i++) {
        int p = primes[i];
        for (int sq: sqrts[i]) {
            // want p | x^2 - n
            // i.e. x^2 \equiv n mod p
            // solutions of x are (s0 + pk), (s1 + pk), where s0 and s1 are square roots of n mod p
            int min_idx = l - (l % p) + sq; // p | (sq + pk)^2 - n for all k >= 0
            if (min_idx < l) min_idx += p;
            assert(min_idx >= l);
            for (int j = min_idx; j <= r; j += p) {
                int idx = j - l;
                assert(idx >= 0 and idx <= r - l);
                while (elts[idx] % p == 0) {
                    elts[idx] /= p;
                    prime_powers[idx][i]++;
                }
            }
        }
    }

    basis.resize(num_primes);
    int cur_basis = 0;
    for (int i = l; i <= r; i++) {
        int idx = i - l;
        if (elts[idx] == 1) {
            // i^2 - n is smooth, add to basis
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
    // cout << "null space size: " << null_space.size() << "\n";
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
                    // cout << "pw: "<< (long long)pw << "\n";
                    assert(pw < max_exp);
                    b_vec[j] += pw;
                }
                a %= n;
            }
        }
        for (int i = 0; i < num_primes; i++) {
            b_vec[i] /= 2;
            int pw = b_vec[i];
            // cout << "pw2: " << (long long)pw << "\n";
            assert(pw < max_exp);
            b *= all_powers[i][pw];
            b %= n;
        }
        // cout << "a: " << (long long)a << ", b: " << (long long)b << "\n";
        // cout << "a^2: " << (long long)((a * a) % n) << ", b^2: " << (long long) ((b * b) % n) << "\n";
        // if((a * a) % n != (b * b) % n) {
        //     cout << "a: " << a << ", b: " << b << "\n";
        //     return;
        // }
        // a^2 = b^2 mod n
        // ll pp = (a + b) % n;
        // ll qq = (a - b) % n;
        // if (qq < 0) qq += n;

        // assert((a * a) % n == (b * b) % n);
        if (a < 0) a += n;
        if (b < 0) b += n;

        ll apb = (a + b) % n;
        ll amb = (a - b) % n;
        if (amb < 0) amb += n;

        // cout << "a: " << a << ", b: " << b << "\n";
        if (abs(a) == abs(b) or a == n - b) continue;
        ll pp = __gcd(n, apb);
        ll qq = __gcd(n, amb);
        if (pp == 1 or qq == 1) continue;


        // cout << "pp: " << pp << ", qq: " << qq << "\n";
        
        // cout << "a: " << a << ", b: " << b << "\n";
        {
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

int2 main() {
    // freopen("A.in", "r", stdin);
    // freopen("A.out", "w", stdout);
    long long nn;
    cin >> nn;
    n = nn;
    m = (ll)sqrt(n);

    // cout << "m: " << m << "\n";

    process_primes();

    find_sqrts(); // find sqrts of n mod p for each p

    sieve(); // find prime factorizations of f(x) for all x in [m + 1, m + sieve_bound]

    find_null_space(); // find solutions to Ax = 0

    try_coeffs(); // iterate over a, b such that a^2 = b^2 mod n

    if (~p) {
        cout << (long long)n << " " << (long long)p << " " << (long long)q << "\n";
    } else {
        cout << "No factorization found!\n";
    }
    
}