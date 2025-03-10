using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.LWE;

public record NTTInfos(
    int n,
    ZnBigInt w,
    KMatrix<ZnBigInt> ntt,
    ZnBigInt[] wPows,
    Rational t,
    KMatrix<ZnBigInt> intt,
    ZnBigInt[] iwPows);