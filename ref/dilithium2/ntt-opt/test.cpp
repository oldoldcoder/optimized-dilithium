#include <vector>
#include <cmath>
#include "iostream"

using namespace std;

const uint32_t N = 256;
uint32_t logN = 8;
uint32_t pVec = 8380417 ;                           //各个素数
uint32_t invp = 58728449  ;                         //pvec的原根的逆原
uint32_t primeLen = 23;                             //pVec的长度
uint32_t inv = 8364049 ;                            //M的逆
uint32_t prVec;

// 模乘运算
void mulMod(uint32_t &r, uint32_t a, uint32_t b, uint32_t m) {
    uint64_t mul = static_cast<uint64_t>(a) * b;
    mul %= static_cast<uint64_t>(m);
    r = static_cast<uint32_t>(mul);
}
// 模指数运算
uint32_t powMod(uint32_t x, uint32_t y, uint32_t modulus) {
    uint32_t res = 1;
    while (y > 0) {
        if (y & 1) {
            mulMod(res, res, x, modulus);
        }
        y = y >> 1;
        mulMod(x, x, x, modulus);
    }
    return res;
}

//素数分解
void get_fac(vector<uint32_t> &s, uint32_t number){
    while (number % 2 == 0) {
        s.push_back(2);
        number /= 2;
    }
    for (uint32_t i = 3; i < sqrt(number); i++) {
        while (number % i == 0) {
            s.push_back(i);
            number /= i;
        }
    }
    if (number > 2) {
        s.push_back(number);
    }
}
// barrett Mul
void mulModBarrett(uint32_t &r, uint32_t a, uint32_t b, uint32_t m, uint32_t pr){
    uint64_t mul = static_cast<uint64_t>(a) * b;
    auto abot = static_cast<uint32_t>(mul);
    auto atop = static_cast<uint32_t>(mul >> 32);
    // int32_t pr = (static_cast<unsigned __int128>(1) << 2*(60+1)) / m;
    uint64_t tmp = static_cast<uint64_t>(abot) * pr;
    tmp >>= 32;
    tmp += static_cast<uint64_t>(atop) * pr;
    tmp >>= 2*(primeLen + 1) - 32;
    tmp *= m;
    tmp = mul - tmp;
    r = static_cast<uint32_t>(tmp);
    if ( r > m ) r -= m;
}

// 模指数运算
uint32_t powModBarrett(uint32_t x, uint32_t y, uint32_t modulus, uint32_t pr) {
    uint32_t res = 1;
    while (y > 0) {
        if (y & 1) {
            mulModBarrett(res, res, x, modulus, pr);
        }
        y = y >> 1;
        mulModBarrett(x, x, x, modulus, pr);
    }
    return res;
}

// 费马小定理求逆元
uint32_t invMod(uint32_t x, uint32_t m) {
    return powMod(x, m - 2, m);
}
// 生成素数原根
uint32_t get_prime_root(uint32_t modulus){
    vector<uint32_t> s;
    uint32_t phi = modulus - 1;
    get_fac(s, phi);
    for (uint32_t r = 2; r <= phi; r++) {
        bool flag = false;
        for (auto it = s.begin(); it != s.end(); it++) {
            if (powModBarrett(r, phi / (*it), modulus, prVec) == 1) {
                flag = true;
                break;
            }
        }
        if (!flag) {
            return r;
        }
    }

    return 0;
}
//  寻找原根
uint32_t findMthRootOfUnity(uint32_t M, uint32_t mod, uint32_t pr) {
    uint32_t res;
    res = get_prime_root(mod);
    if((mod - 1) % M == 0) {
        uint32_t factor = (mod - 1) / M;
        res = powModBarrett(res, factor, mod, pr);
        return res;
    }else{
        return 0;
    }
}
uint32_t REDC(uint32_t a, uint32_t b){
    uint64_t U = static_cast<uint64_t>(a) * b;
    auto U0 = static_cast<uint64_t>(U);
    uint32_t U1 = U >> 32;
    uint32_t Q = U0 * invp;
    uint64_t Hx = static_cast<uint64_t>(Q) * pVec;
    uint32_t H = Hx >> 32;
    uint32_t V = U1 < H ? U1 + pVec - H : U1 - H;
    return V;
}


int generate_prime(){
    uint32_t* g;                               //pvec的原根
    uint32_t* gi;                              //pvec的原根的逆原

    // 这里这个除法除不干净
    prVec = static_cast<::uint32_t>((static_cast<::uint64_t>(1) << 2 * (primeLen + 1)) / pVec);

    g = new uint32_t[N >> 1];
    gi = new uint32_t[N >> 1];
    cout << "prVec:" << prVec << endl;
    uint32_t root = findMthRootOfUnity(N, pVec, prVec);   // get prime N-root of mod
    cout << "root:" << root << endl;
    uint32_t invroot = invMod(root, pVec);                // get inv prime N-root of mod
    cout << "invroot:" << invroot << endl;

    uint32_t temp = static_cast<uint32_t>(1) << 16;
    g[0] = 1;
    mulModBarrett(g[0], g[0], temp, pVec, prVec);
    mulModBarrett(g[0], g[0], temp, pVec, prVec);
    gi[0] = g[0];

    for (int P = 1; P < N >> 1; P++) {
        g[P] = REDC(g[P - 1], root);
        mulModBarrett(g[P], g[P], temp, pVec, prVec);
        mulModBarrett(g[P], g[P], temp, pVec, prVec);
        gi[P] = REDC(gi[P - 1], invroot);
        mulModBarrett(gi[P], gi[P], temp, pVec, prVec);
        mulModBarrett(gi[P], gi[P], temp, pVec, prVec);
    }

    mulModBarrett(inv, inv, temp, pVec, prVec);
    mulModBarrett(inv, inv, temp, pVec, prVec);
    // 打印具体的g和gi的值
    cout<<"g table :";
    for(int i = 0 ; i < N / 2 ; ++i){
        cout<< g[i] << ",";
    }
    cout << endl;
    cout<<"gi table :";
    for(int i = 0 ; i < N / 2 ; ++i){
        cout<< gi[i] << ",";
    }

    return 1;
}

int rev[N];

// bit反转，蝶形变换
void bitreverse(){
    int revint[N];
    for (int i = 0; i < N; ++i)
    {
        revint[i] = 0;
        rev[i] = i;
    }
    for (int i = 0; i < N; ++i){				// get rev
        revint[i] = (revint[i >> 1] >> 1) | ((i & 1) << (logN - 1));
        if(i < revint[i]){
            rev[revint[i]] = i;
            rev[i] = revint[i];
        }
    }
}
int main(){
    generate_prime();
}