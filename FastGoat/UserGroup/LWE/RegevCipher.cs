using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.LWE;

public readonly struct RegevCipher
{
    public Vec<ZnBigInt> A { get; }
    public ZnBigInt B { get; }

    public RegevCipher(Vec<ZnBigInt> a, ZnBigInt b)
    {
        (A, B) = (a, b);
    }

    public override string ToString() => $"[{A}, {B}]";

    public static implicit operator RegevCipher((Vec<ZnBigInt> a, ZnBigInt b) e) => new(e.a, e.b);
    public static RegevCipher operator +(RegevCipher a, RegevCipher b) => new(a.A + b.A, a.B + b.B);
    public static RegevCipher operator +(RegevCipher a, int b) => new(a.A, a.B + b);
}