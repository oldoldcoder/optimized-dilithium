#include <iostream>
#include <vector>
#include <cstdint>

using namespace std;

uint32_t mod(uint32_t a, uint32_t p) {
    return (a % p + p) % p;
}

void validateRoots(uint32_t root, uint32_t invroot, uint32_t pVec, int N) {
    // 验证原根
    vector<uint32_t> powers(N);
    powers[0] = 1; // g^0
    for (int i = 1; i < N; ++i) {
        powers[i] = mod(powers[i - 1] * root, pVec);
    }

    cout << "Checking original root powers:" << endl;
    for (int i = 0; i < N; ++i) {
        cout << powers[i] << " ";
    }
    cout << endl;

    // 验证逆元根
    uint32_t product = mod(root * invroot, pVec);
    cout << "Checking if root * invroot == 1: " << product << endl;

    if (product == 1) {
        cout << "Validation successful: root and invroot are correct." << endl;
    } else {
        cout << "Validation failed: root and invroot are incorrect." << endl;
    }
}

int main() {
    uint32_t root = 6644104/* 从你的代码中获得原根 */;
    uint32_t invroot = 6125690/* 从你的代码中获得逆原根 */;
    uint32_t pVec = 8380417/* 模数 */;
    int N = 256/* N的值 */;

    validateRoots(root, invroot, pVec, N);
    return 0;
}
