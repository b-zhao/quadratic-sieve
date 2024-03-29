# quadratic-sieve

## How to run the quadratic sieve
* The program reads n from stdin and outputs its two prime factors p and q
* Complilation:
    g++ qs.cpp -o qs, Running: ./qs (then input n)

## Stress tester
* Stress tester taken directly from https://github.com/FlowerOfSorrow/Utility/blob/main/compare.cpp
* Compile gen.cpp, qs.cpp, naive.cpp, compare.cpp and run ./compare naive qs

## Progress

* Passes the stress tester for most n = pq where p and q are odd (picked randomly) and p, q < 1e6
* TODO: Determine why there seems to be no solution for some n = pq, where max(p, q) > 1e6

Sources used:
https://www.cc.gatech.edu/~rpeng/CS6550_S21/Feb15QS.pdf, https://www.cs.virginia.edu/crab/QFS_Simple.pdf, https://en.wikipedia.org/wiki/Quadratic_sieve
