// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stack>

using namespace std;
using namespace Rcpp;

double euclidean_distance(NumericVector x, NumericVector y) {
    return sqrt(sum(pow(x - y, 2)));
}

double manhattan_distance(NumericVector x, NumericVector y) {
    return sum(abs(x - y));
}

typedef double (*distance_function_pointer)(NumericVector, NumericVector);

distance_function_pointer select_distance(string distance) {
    distance_function_pointer dist_ptr;
    if (distance == "euclidean")
        dist_ptr = &euclidean_distance;
    else if (distance == "manhattan")
        dist_ptr = &manhattan_distance;
    else
        stop("Unknown distance!");
    return dist_ptr;
}

NumericMatrix distances(NumericMatrix X, std::string distance) {
    distance_function_pointer dist_ptr = select_distance(distance);
    int n = X.nrow();
    NumericMatrix D(n, n);
    for (int i = 0; i < n; i++)
        for(int j = i + 1; j < n; j++) {
            D(i, j) = (*dist_ptr)(X.row(i), X.row(j));
            D(j, i) = D(i, j);
        }
    return D;
}

IntegerVector order(NumericVector x) {
    return match(clone(x).sort(), x) - 1;
}

IntegerMatrix neighbors_among_training(NumericMatrix D) {
    // row = idx of sample unit; col (r - 1) = idx of r-th neighbor
    IntegerMatrix N(D.nrow(), D.ncol() - 1);
    for (int i = 0; i < D.nrow(); i++)
        N(i, _) = order(D.row(i))[Range(1, D.ncol() - 1)];
    return N;
}

NumericMatrix test_training_distances(NumericMatrix X_tst, NumericMatrix X_trn, std::string distance) {
    distance_function_pointer dist_ptr = select_distance(distance);
    int n_tst = X_tst.nrow();
    int n_trn = X_trn.nrow();
    NumericMatrix D(n_tst, n_trn);
    for (int i = 0; i < n_tst; i++)
        for (int j = 0; j < n_trn; j++)
            D(i, j) = (*dist_ptr)(X_tst.row(i), X_trn.row(j));
    return D;
}

IntegerMatrix neighbors(NumericMatrix D) {
    IntegerMatrix N(D.nrow(), D.ncol() - 1);
    for (int i = 0; i < D.nrow(); i++) N(i, _) = order(D.row(i));
    return N;
}

IntegerMatrix neighbor_of(NumericMatrix D, NumericMatrix D_trn, IntegerMatrix N_trn) {
    int n_tst = D.nrow();
    int n_trn = D_trn.nrow();
    IntegerMatrix N_of(n_trn, n_tst);
    // row = trn; col = tst
    // { cell value = r } <=> { tst is the r-th neighbor of trn }
    for (int i = 0; i < n_trn; i++)
        for (int j = 0; j < n_tst; j++)
            for (int r = 1; r <= n_trn - 1; r++)
                if (D(j, i) <= D_trn(i, N_trn(i, r - 1))) {
                    N_of(i, j) = r;
                    break;
                }
    return N_of;
}

IntegerMatrix count_simple_cycles(IntegerMatrix N) {
    int n = N.nrow();
    // row: r = 1, ... , n - 1; col: m = 2, ... , n
    IntegerMatrix C(n - 1, n);
    for (int r = 1; r <= n - 1; r++) {
        map<int, bool> visited;
        for (int i = 0; i < n; i++) { // O(n)
            if (visited[i]) continue;
            int j = i;
            stack<int> walk;
            while (!visited[j]) {
                visited[j] = true;
                walk.push(j);
                j = N(j, r - 1);
            }
            int cycle_size = 1;
            while (!walk.empty()) {
                cycle_size++;
                walk.pop();
                if ((!walk.empty()) && (walk.top() == j)) {
                    C(r - 1, cycle_size - 1)++;
                    break;
                }
            }
        }
        // the following is used for the "first term" in Z_r(\beta_r)
        C(r - 1, 0) = n; 
        for (int m = 2; m <= n; m++)
            C(r - 1, 0) -= m * C(r - 1, m - 1);
    }
    return C;
}

IntegerVector same_class_neighbors(IntegerVector y, IntegerMatrix N) {
    int n = y.length();
    IntegerVector A(n - 1);
    for (int r = 1; r <= n - 1; r++)
        for (int i = 0; i < n; i++)
            if (y[i] == y[N(i, r - 1)]) A[r - 1]++;
    return A;
}

double log_Z_r(double beta_r, IntegerVector C_r, int L) {
    int n = C_r.length();
    double res = C_r[0] * (beta_r + log(exp(-beta_r) * (L - 1) + 1));
    for (int m = 2; m <= n; m++) {
         if (C_r[m - 1] == 0) continue;
         res += C_r[m - 1] * (m*beta_r + log(pow(exp(-beta_r) * (L - 1) + 1, m) 
                                             + (L - 1) * pow(1 - exp(-beta_r), m)));
    }
    return res;
}

// [[Rcpp::export(.log_p_r)]]
double log_p_r(double beta_r, int A_r, IntegerVector C_r, int L) {
    return beta_r * A_r - log_Z_r(beta_r, C_r, L);
}

// [[Rcpp::export(.data_topology)]]
List data_topology(NumericMatrix X_trn, IntegerVector y_trn, std::string distance) {
    List dt;

    NumericMatrix D = distances(X_trn, distance);
    IntegerMatrix N = neighbors_among_training(D);
    IntegerMatrix C = count_simple_cycles(N);
    IntegerVector A = same_class_neighbors(y_trn, N);
    
    dt["X_trn"] = X_trn;
    dt["y_trn"] = y_trn;
    dt["distance"] = distance;
    dt["L"] = unique(y_trn).length();
    dt["D_trn"] = D;
    dt["N_trn"] = N;
    dt["C"] = C;
    dt["A"] = A;

    return dt;
}

// e^{q_i} / \sum_j e^{q_j} = exp(q_i - M - log(\sum_j e^{q_j - M})), where M = \max_i q_i

// [[Rcpp::export(.loocv)]]
List loocv(List dt, NumericVector beta_hat) {
    List cv;
    NumericVector y = dt["y_trn"];
    int n = y.length();
    int L = dt["L"];
    IntegerMatrix N = dt["N_trn"];
    int k_max = beta_hat.length();
    
    arma::cube P(k_max, n, L); // , arma::fill::zeros);
    for (int r = 1; r <= k_max; r++) {
        for (int i = 0; i < n; i++) {
            NumericVector q(L);
            for (int ell = 0; ell < L; ell++) {
                if (ell == y[N(i, r - 1)]) q[ell]++;
                for (int j = 0; j < n; j++) {
                     if (j == i) continue;
                     if (N(j, r - 1) == i && y[j] == ell) q[ell]++;
                }
            }
            q = beta_hat[r - 1] * q;
            double M = max(q);
            double log_sum_exp = log(sum(exp(q - M)));
            for (int ell = 0; ell < L; ell++)
                P(r - 1, i, ell) = exp(q[ell] - M - log_sum_exp);
        }
    }

    NumericMatrix Q(n, n - 1);
    for (int r = 1; r <= k_max; r++)
        for (int i = 0; i < n; i++)
            Q(i, r - 1) = P(r - 1, i, y[i]);
    cv["Q"] = Q;
    
    NumericVector S(k_max); // accuracy
    for (int k = 1; k <= k_max; k++) {
         for (int i = 0; i < n; i++) {
             NumericVector pr(L);
             for (int r = 1; r <= k; r++)
                 for (int ell = 0; ell < L; ell++)
                     pr[ell] += P(r - 1, i, ell);
             pr = pr / k;
             if (which_max(pr) == y[i]) S[k - 1]++;
         }
         S[k - 1] /= n;
    }
    cv["error"] = 1.0 - double(max(S));
    cv["k_hat"] = which_max(S) + 1;
    
    return cv;
}

// [[Rcpp::export(.predict_test)]]
List predict_test(List model, NumericMatrix X_tst) {
    List pred;
    int n_tst = X_tst.nrow();
    List dt = model["dt"];
    NumericVector y_trn = dt["y_trn"];
    int L = dt["L"];
    NumericVector beta_hat = model["beta_hat"];
    int k_hat = model["k_hat"];

    NumericMatrix D = test_training_distances(X_tst, dt["X_trn"], dt["distance"]);
    IntegerMatrix N = neighbors(D);
    IntegerMatrix N_of = neighbor_of(D, dt["D_trn"], dt["N_trn"]);
    int n_trn = N_of.nrow();

    NumericMatrix prob(n_tst, L);
    for (int i = 0; i < n_tst; i++) {
        for (int r = 1; r <= k_hat; r++) {
            NumericVector q(L);
            for (int ell = 0; ell < L; ell++) {
                if (ell == y_trn[N(i, r - 1)]) q[ell]++;
                for (int j = 0; j < n_trn; j++)
                    if (N_of(j, i) == r && y_trn[j] == ell) q[ell]++;
            }
            q = beta_hat[r - 1] * q;
            double M = max(q);
            double log_sum_exp = log(sum(exp(q - M)));
            prob(i, _) = prob(i, _) + exp(q - M - log_sum_exp);
        }
        prob(i, _) = prob(i, _) / k_hat;
    }
    
    pred["X_tst"] = X_tst;
    pred["D"] = D;
    pred["N"] = N;
    pred["N_of"] = N_of;
    pred["prob"] = prob;

    return pred;
}
