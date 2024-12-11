using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Lattice;

public struct CypherLWE(ZnInt64[] ai, ZnInt64 b)
{
    public ZnInt64[] Ai => ai;
    public ZnInt64 B => b;

    public static CypherLWE Not(CypherLWE A)
    {
        var Ai = A.Ai.Select(e => e.Opp()).ToArray();
        var Bi = A.B.Opp();
        return new(Ai, Bi);
    }

    public static CypherLWE Xor(CypherLWE A, CypherLWE B)
    {
        var Ai = A.Ai.Zip(B.Ai).Select(e => e.First + e.Second).ToArray();
        var Bi = A.B + B.B;
        return new(Ai, Bi);
    }
}

public struct CypherRLWE(EPoly<ZnInt> a, EPoly<ZnInt> b)
{
    public EPoly<ZnInt> A => a;
    public EPoly<ZnInt> B => b;

    public static CypherRLWE Not(CypherRLWE C) => new(-C.A, -C.B);

    public static CypherRLWE Xor(CypherRLWE C0, CypherRLWE C1) => new(C0.A + C1.A, C0.B + C1.B);
}