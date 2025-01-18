using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.LWE;

public readonly struct RegevCipher
{
    public Vec<ZnInt64> A { get; }
    public ZnInt64 B { get; }

    public RegevCipher(Vec<ZnInt64> a, ZnInt64 b)
    {
        (A, B) = (a, b);
    }

    public override string ToString() => $"[{A}, {B}]";

    public static implicit operator RegevCipher((Vec<ZnInt64> a, ZnInt64 b) e) => new(e.a, e.b);
    public static RegevCipher operator +(RegevCipher a, RegevCipher b) => new(a.A + b.A, a.B + b.B);
    public static RegevCipher operator +(RegevCipher a, int b) => new(a.A, a.B + b);
}