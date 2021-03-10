# quadratic-sieve

## How to run the quadratic sieve
* The program reads n from stdin and outputs its two prime factors p and q
* Complilation:
    g++ qs.cpp -o qs, Running: ./qs (then input n)

## Stress tester
* Stress tester taken directly from https://github.com/FlowerOfSorrow/Utility/blob/main/compare.cpp
* Compile gen.cpp, qs.cpp, naive.cpp, compare.cpp and run ./compare naive qs

## Progress

* Passes the stress tester for most n = pq where p and q are odd (picked randomly) and p, q < 10000
* Some cases where the factorization is not found: n = 5026163 = 601 * 8363
* TODO: Determine a better way to iterate through / sample vectors from the null space? (i.e. the space that gives a, b such
that a^2 = b^2 mod n)

Sources used:
https://www.cc.gatech.edu/~rpeng/CS6550_S21/Feb15QS.pdf, https://www.cs.virginia.edu/crab/QFS_Simple.pdf, https://en.wikipedia.org/wiki/Quadratic_sieve
