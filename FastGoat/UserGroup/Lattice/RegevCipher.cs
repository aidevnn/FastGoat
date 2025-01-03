using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Lattice;

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
}