#define sparse_binop(_name, _op) \
int _name (sil_State *S) { \
    const SparseMat *a = (const SparseMat *)sil_topointer(S, 1); \
    const SparseMat *b = (const SparseMat *)sil_topointer(S, 2); \
    if(a == NULL || b == NULL) \
        return sil_err(S, "Invalid arguments"); \
    SparseMat *c = new SparseMat(*a _op *b); \
    sil_settop(S, 0); \
    sil_pushSparse(S, c); \
    return 0; \
}

// create add, mul, and sub functions
sparse_binop(add, +);
sparse_binop(sub, -);

SparseMat mul_slow(sil_State *S, int t1, int t2) {
    if(t1 == 5 || t1 == 6 || t1 == 7) {
        double x = sil_todouble(S, 1);
        const SparseMat *y = (const SparseMat *)sil_topointer(S, 2);
        return x * (*y);
    }
    const SparseMat *x = (const SparseMat *)sil_topointer(S, 2);
    double y = sil_todouble(S, 1);
    return (*x) * y;
}

extern "C" {

int mul(sil_State *S) {
    int t1 = sil_type(S, 1);
    int t2 = sil_type(S, 2);
    if(t1 != 16 || t2 != 16) { /* slow version */
        SparseMat *out = new SparseMat(mul_slow(S, t1, t2));
        sil_settop(S, 0);
        sil_pushSparse(S, out);
        return 0;
    }
    const SparseMat *a = (const SparseMat *)sil_topointer(S, 1);
    const SparseMat *b = (const SparseMat *)sil_topointer(S, 2);
    if(a == NULL || b == NULL)
        return sil_err(S, "Invalid arguments");
    SparseMat *c = new SparseMat(a->cwiseProduct(*b));
    sil_settop(S, 0);
    sil_pushSparse(S, c);
    return 0;
}

// simple matrix creation routines
int eye(sil_State *S) {
    if(sil_type(S, 1) != 6) return sil_err(S, "eye requires int argument");
    int n = sil_tointeger(S, 1);
    if(n < 0 || n > 100000) {
        return sil_err(S, "eye invalid size");
    }
    SparseMat *I = new SparseMat(n,n);

    typedef Triplet<double> T;
    std::vector<T> trip;
    trip.reserve(n);

    for(int j=0; j<n; j++) {
        trip.push_back(T(j,j,1.0));
    }

    I->setFromTriplets(trip.begin(), trip.end());
    sil_settop(S, 0);
    sil_pushSparse(S, I);
    return 0;
}

// [(Int,Int,Float)] -> SparseMat
int fromList(sil_State *S) {
    int k = sil_llen(S, 1);
    int n = sil_tointeger(S, 2);
    int m = sil_tointeger(S, 3);
    typedef Triplet<double> T;
    std::vector<T> trip;

    for(int e=0; e<k; e++) {
        sil_settop(S, 1); // 1 is top elem
        sil_behead(S, 1);
        if(sil_unpack(S, 2) != 3) {
            return sil_err(S, "Improperly formatted list.\n");
        }
        trip.push_back(T(sil_tointeger(S, 3),
                         sil_tointeger(S, 4),
                         sil_todouble(S, 5)));
    }
    SparseMat *A = new SparseMat(n,m);

    A->setFromTriplets(trip.begin(), trip.end());
    sil_settop(S, 0);
    sil_pushSparse(S, A);
    return 0;
}

// SparseMat v -> [(Int,Int,Float)]
int toList(sil_State *S) {
    const SparseMat *A = (const SparseMat *)sil_topointer(S, 1);
    int n = 0;
    if(A == NULL) {
        return sil_err(S, "toList: requires SparseMat");
    }
    for(int j=0; j<A->outerSize(); ++j) {
        for(SparseMat::InnerIterator it(*A,j); it; ++it) {
            sil_pushinteger(S, it.row());
            sil_pushinteger(S, j);
            sil_pushdouble(S, it.value());
            sil_settuple(S, 3);
            n++;
        }
    }
    sil_setcons(S, n);
    sil_remove(S, 1); // clear mat
    return 0;
}

#define READ_WRAP(out, pos, name, max) \
    int out = sil_tointeger(S, pos); \
    if(out < 0) out += max; \
    if(out < 0 || out >= max) \
        return sil_err(S, "Invalid " #name " %d", out); \

// SparseMat v, Int i, Int j -> Float
int elem(sil_State *S) {
    const SparseMat *A = (const SparseMat *)sil_topointer(S, 1);
    if(A == NULL) {
        return sil_err(S, "Invalid argument");
    }

    READ_WRAP(i, 2, "row", A->rows());
    READ_WRAP(j, 3, "column", A->cols());
    sil_settop(S, 1);
    // FIXME: implement own binary search / return whole col?
    double x = ((SparseMat *)A)->coeffRef(i,j);
    sil_pushdouble(S, x);
    sil_remove(S, 1);
    return 0;
}

// Int i, Int j, Float x -> ST(SparseMat, Nil)
int setElem(sil_State *S) {
    size_t len; // Always assume len is wrong!
    SparseMat *A = (SparseMat *)sil_getST(S, &len);
    if(A == NULL) {
        return sil_err(S, "Can't update - no SparseMat present");
    }
    READ_WRAP(i, 2, "row", A->rows());
    READ_WRAP(j, 3, "column", A->cols());
    A->coeffRef(i,j) = sil_todouble(S, 4);

    sil_settop(S, 0);
    sil_pushnil(S);
    return 0;
}

// (Int col -> [(Int,Int,Float)]), Int n, Int m -> SparseMat
int fromColFn(sil_State *S) {
    int n = sil_tointeger(S, 2);
    int m = sil_tointeger(S, 3);
    if(n < 0 || n > 10000000 || m < 0 || m > 10000000) {
        return sil_err(S, "invalid matrix size (%d, %d)", n, m);
    }
    typedef Triplet<double> T;
    std::vector<T> trip;

    for(int j=0; j<m; j++) {
        sil_settop(S, 1);
        sil_pushvalue(S, 1); // copy function
        sil_pushinteger(S, j); // push arg
        sil_call(S, 2);

        int entries = sil_unpack(S, 2);
        for(int k=0; k<entries; k++) {
            if(sil_unpack(S, 2) != 2) {
                return sil_err(S, "invalid function return value");
            }
            trip.push_back(T(sil_tointeger(S, 2), j, sil_todouble(S, 3)));
            sil_remove(S, 2); sil_remove(S, 2); sil_remove(S, 2);
        }
    }
    sil_settop(S, 0);

    SparseMat *A = new SparseMat(n,m);
    A->setFromTriplets(trip.begin(), trip.end());
    sil_pushSparse(S, A);
    return 0;
}

}
