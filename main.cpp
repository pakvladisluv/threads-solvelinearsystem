#include "func.h"

typedef struct _ARGS {
    double *matrix = NULL;
    double *vector = NULL;
    double *result = NULL;
    double *tmp_common = NULL;
    double *trash = NULL;
    int n;
    int thread_num;
    int total_threads;
} ARGS;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
bool check = true;

void *Cholesky(void *p_arg) {
    ARGS *arg = (ARGS *)p_arg;
    unsigned long long time;
    time = currentTimeNano();
    solving(arg->n, arg->matrix, arg->tmp_common, arg->vector, arg->result, arg->trash, arg->thread_num, arg->total_threads);
    pthread_mutex_lock(&mutex);
    time = currentTimeNano() - time;
    cout << "Time of thread " << arg->thread_num << ": " << time << " ns" << endl;
    pthread_mutex_unlock(&mutex);
    return NULL;
}

int main(int argc, char *argv[]) {
    if (argc < FIVE || argc > SIX) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    int p = atoi(argv[1]);
    int n = atoi(argv[2]);
    int m = atoi(argv[3]);
    int k = atoi(argv[4]);
    char *name = argv[FIVE];
    if (k == 0 && argc <= FIVE) {
        return -1;
    }
    if (n < 1 || m < 1 || k < 0 || k > 4 || n < m) {
        return -1;
    }
    if (p < 1) {
        return -1;
    }
    if (n > MAX) {
        return -2;
    }
    double *a = new double[n * n];
    double *R = new double[n * n];
    if (k == 0) {
        if (!read(name, a, n)) {
            delete[] a;
            delete[] R;
            return -3;
        }
    } else {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                a[i * n + j] = func(k, n, i, j);
                R[i * n + j] = 0;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (fabs(a[i * n + j] - a[j * n + i]) > MIN) {
                delete[] a;
                delete[] R;
                return -4;
            }
        }
    }
    double *b = new double[n];
    right(a, b, n);
    print(n, m, a, b);
    cout << endl;
    if (k == 3) {
        delete[] a;
        delete[] R;
        delete[] b;
        return 0;
    }
    for (int i = 0; i < n; i++) {
        if(fabs(a[i * n + i]) < MIN){
            delete[] a;
            delete[] R;
            delete[] b;
            return -4;
        }
    }
    double *D = new double[n];
    double *X = new double[n];
    for (int i = 0; i < n; i++) {
        D[i] = 0;
    }
    ARGS *args;
    args = new ARGS[p];
    pthread_t *threads;
    threads = new pthread_t[p];
    for (int i = 0; i < p; i++) {
        args[i].matrix = a;
        args[i].vector = b;
        args[i].result = D;
        args[i].tmp_common = R;
        args[i].trash = X;
        args[i].n = n;
        args[i].thread_num = i;
        args[i].total_threads = p;
    }
    unsigned long long time = 0;
    time = currentTimeNano();
    for (int i = 0; i < p; i++) {
        if (pthread_create(threads + i, 0, Cholesky, args + i) != 0) {
            cout << "Can not create thread " << i << '!' << endl;
            delete[] a;
            delete[] b;
            delete[] R;
            delete[] D;
            delete[] X;
            delete[] threads;
            delete[] args;
            return -FIVE;
        }
    }
    for (int i = 0; i < p; i++) {
        if (pthread_join(threads[i], 0) != 0) {
            cout << "Can not wait thread " << i << '!' << endl;
            delete[] a;
            delete[] b;
            delete[] R;
            delete[] D;
            delete[] X;
            delete[] threads;
            delete[] args;
            return -FIVE;
        }
    }
    delete[] threads;
    time = currentTimeNano() - time;
    cout << "Time: " << time << endl;
    if (!determinant(args->n, args->tmp_common)) {
        delete[] a;
        delete[] b;
        delete[] R;
        delete[] D;
        delete[] X;
        delete[] args;
        return -4;
    }
    cout << "Solution:";
    for (int i = 0; i < m; i++) {
        cout << ' ' << args->result[i];
    }
    cout << endl;
    residual(args->n, args->matrix, args->vector, args->result);
    error(n, D);
    delete[] a;
    delete[] b;
    delete[] R;
    delete[] D;
    delete[] X;
    delete[] args;
    return 0;
}
