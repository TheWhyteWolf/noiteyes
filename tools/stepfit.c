/* stepfit.c -- anneal a hidden order pi for the short-step gak41 model.
 *
 * Model: c_t = D[a_t], a_t = a_(t-1) + b_(t-1) * k(p_t) mod 83, steps k in
 * 1..K, b in the order-41 subgroup (quadratic residues), piecewise constant
 * (it changes only at twist letters). With pi = D^-1 every consecutive
 * difference d_t = pi(c_t) - pi(c_(t-1)) must lie in b_t * [1, K].
 *
 * Cost(pi) = sum over messages of a Viterbi pass over b (41 states):
 *   1 per position whose d_t is not in b*[1,K], mu per change of b,
 *   mu if the first b is not 1 (the state starts at b = 1).
 * Positions 0 (off-chain header) and 1 (difference to the header) are skipped.
 * Simulated annealing over pi by transpositions; best pi printed.
 *
 * Twist-tied mode (NT > 0): b may change only at a position whose step
 * equals one of NT twist steps, by that twist's multiplier (cost_tied);
 * other b jumps cost mu. Twist (step, mult) pairs are given or annealed.
 *
 * stdin: K mu iters seed T0 T1 / NT (step mult)*NT [0 0 = anneal] /
 *        nmsg / (len c_0 .. c_len-1) per message
 *        [iters == 0: then 83 ints pi, which is only scored]
 * stdout: best_cost / 83 ints pi
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define P 83
#define NB 41
#define MAXM 16
#define MAXL 256

static int nm, len[MAXM], msg[MAXM][MAXL];
static int QRv[NB], inb[P][NB];
static int K, mu, NT;              /* NT > 0: twist-tied mode */
static int tstep[8], tmul[8];      /* twist letter j: step value, multiplier index */
static int qidx[P], mulidx[NB][NB]; /* index of a QR; index of QRv[i]*QRv[j] */

static unsigned long long rs;
static inline unsigned long long rnd(void) {
    rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs;
}
static inline double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }

/* Twist-tied cost: b changes ONLY by multiplier QRv[tmul[j]] at a position
 * whose step s = d * b^-1 equals tstep[j] (the twist letter's own step). A
 * position with s outside [1,K] costs 1; any other b jump costs mu. */
static int cost_tied(const int *pi) {
    int total = 0, prev[NB], cur[NB], nxt[NB], c[NB];
    static int sof[P][NB], init = 0, invq[NB];
    if (!init) {
        for (int j = 0; j < NB; j++) {
            int inv = 1;
            for (int e = 0; e < P - 2; e++) inv = inv * QRv[j] % P;
            invq[j] = inv;
        }
        for (int d = 0; d < P; d++) for (int j = 0; j < NB; j++) sof[d][j] = d * invq[j] % P;
        init = 1;
    }
    for (int m = 0; m < nm; m++) {
        int L = len[m];
        if (L < 3) continue;
        for (int j = 0; j < NB; j++) prev[j] = j == 0 ? 0 : mu;
        for (int t = 2; t < L; t++) {
            int d = ((pi[msg[m][t]] - pi[msg[m][t - 1]]) % P + P) % P;
            int best = 1 << 29;
            for (int j = 0; j < NB; j++) {
                int s = sof[d][j];
                c[j] = prev[j] + !(s >= 1 && s <= K);
                nxt[j] = j;
                for (int q = 0; q < NT; q++) if (s == tstep[q]) { nxt[j] = mulidx[j][tmul[q]]; break; }
                if (c[j] < best) best = c[j];
            }
            for (int j = 0; j < NB; j++) cur[j] = best + mu;
            for (int j = 0; j < NB; j++) if (c[j] < cur[nxt[j]]) cur[nxt[j]] = c[j];
            memcpy(prev, cur, sizeof prev);
        }
        int best = prev[0];
        for (int j = 1; j < NB; j++) if (prev[j] < best) best = prev[j];
        total += best;
    }
    return total;
}

static int cost_free(const int *pi) {
    int total = 0, prev[NB], cur[NB];
    for (int m = 0; m < nm; m++) {
        int L = len[m];
        if (L < 3) continue;
        for (int t = 2; t < L; t++) {
            int d = ((pi[msg[m][t]] - pi[msg[m][t - 1]]) % P + P) % P;
            if (t == 2) {
                for (int j = 0; j < NB; j++) prev[j] = (!inb[d][j]) + (j == 0 ? 0 : mu);
                continue;
            }
            int best = prev[0];
            for (int j = 1; j < NB; j++) if (prev[j] < best) best = prev[j];
            best += mu;
            for (int j = 0; j < NB; j++) {
                int s = prev[j] < best ? prev[j] : best;
                cur[j] = s + !inb[d][j];
            }
            memcpy(prev, cur, sizeof prev);
        }
        int best = prev[0];
        for (int j = 1; j < NB; j++) if (prev[j] < best) best = prev[j];
        total += best;
    }
    return total;
}

static int cost(const int *pi) { return NT ? cost_tied(pi) : cost_free(pi); }

int main(void) {
    long iters; unsigned long long seed; double T0, T1;
    if (scanf("%d %d %ld %llu %lf %lf", &K, &mu, &iters, &seed, &T0, &T1) != 6) return 1;
    if (scanf("%d", &NT) != 1 || NT > 8) return 1;
    int tin[8][2];
    for (int q = 0; q < NT; q++) if (scanf("%d %d", &tin[q][0], &tin[q][1]) != 2) return 1;
    if (scanf("%d", &nm) != 1 || nm > MAXM) return 1;
    for (int m = 0; m < nm; m++) {
        if (scanf("%d", &len[m]) != 1 || len[m] > MAXL) return 1;
        for (int i = 0; i < len[m]; i++) if (scanf("%d", &msg[m][i]) != 1) return 1;
    }
    /* QRv[0] = 1 (the start state); inb[d][j]: d * QRv[j]^-1 in [1,K] */
    int n = 0, seen[P] = {0};
    QRv[n++] = 1; seen[1] = 1;
    for (int x = 2; x < P; x++) { int q = x * x % P; if (!seen[q]) { seen[q] = 1; QRv[n++] = q; } }
    for (int j = 0; j < NB; j++) {
        int inv = 1;
        for (int e = 0; e < P - 2; e++) inv = inv * QRv[j] % P;
        for (int d = 0; d < P; d++) { int s = d * inv % P; inb[d][j] = (s >= 1 && s <= K); }
    }
    for (int j = 0; j < NB; j++) qidx[QRv[j]] = j;
    for (int i = 0; i < NB; i++) for (int j = 0; j < NB; j++) mulidx[i][j] = qidx[QRv[i] * QRv[j] % P];
    rs = seed * 2654435761ULL + 88172645463325252ULL;
    for (int q = 0; q < NT; q++) {          /* given (step, mult) or random if 0 */
        tstep[q] = tin[q][0] ? tin[q][0] : 1 + (int)(rnd() % K);
        tmul[q] = tin[q][1] ? qidx[tin[q][1]] : 1 + (int)(rnd() % (NB - 1));
    }
    int pi[P], bestpi[P];
    if (iters == 0) {
        for (int i = 0; i < P; i++) if (scanf("%d", &pi[i]) != 1) return 1;
        printf("%d\n", cost(pi));
        return 0;
    }
    for (int i = 0; i < P; i++) pi[i] = i;
    for (int i = P - 1; i > 0; i--) { int j = rnd() % (i + 1); int t = pi[i]; pi[i] = pi[j]; pi[j] = t; }
    int c = cost(pi), bestc = c, bst[8], bsm[8];
    int fixtw = NT && tin[0][0] != 0;        /* twist params given: keep them fixed */
    memcpy(bestpi, pi, sizeof pi); memcpy(bst, tstep, sizeof bst); memcpy(bsm, tmul, sizeof bsm);
    double lr = log(T1 / T0);
    for (long it = 0; it < iters; it++) {
        double T = T0 * exp(lr * it / iters);
        if (NT && fixtw == 0 && rnd() % 20 == 0) {   /* twist-parameter move */
            int q = rnd() % NT, os = tstep[q], om = tmul[q];
            if (rnd() & 1) tstep[q] = 1 + (int)(rnd() % K); else tmul[q] = 1 + (int)(rnd() % (NB - 1));
            int c2 = cost(pi);
            if (c2 <= c || urand() < exp((c - c2) / T)) {
                c = c2;
                if (c < bestc) { bestc = c; memcpy(bestpi, pi, sizeof pi); memcpy(bst, tstep, sizeof bst); memcpy(bsm, tmul, sizeof bsm); }
            } else { tstep[q] = os; tmul[q] = om; }
            continue;
        }
        int x = rnd() % P, y = rnd() % P;
        if (x == y) continue;
        int t = pi[x]; pi[x] = pi[y]; pi[y] = t;
        int c2 = cost(pi);
        if (c2 <= c || urand() < exp((c - c2) / T)) {
            c = c2;
            if (c < bestc) { bestc = c; memcpy(bestpi, pi, sizeof pi); memcpy(bst, tstep, sizeof bst); memcpy(bsm, tmul, sizeof bsm); }
        } else { t = pi[x]; pi[x] = pi[y]; pi[y] = t; }
    }
    printf("%d\n", bestc);
    for (int i = 0; i < P; i++) printf("%d ", bestpi[i]);
    printf("\n");
    for (int q = 0; q < NT; q++) printf("%d %d ", bst[q], QRv[bsm[q]]);
    printf("\n");
    return 0;
}
