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

public struct CypherRLWE
{
    public CypherRLWE(EPoly<ZnInt> a0, EPoly<ZnInt> p0)
    {
        (a, p) = (a0, p0);
    }
    public EPoly<ZnInt> a { get; }
    public EPoly<ZnInt> p { get; }

    public void Deconstruct(out EPoly<ZnInt> a0, out EPoly<ZnInt> p0) => (a0, p0) = (a, p);

    public static CypherRLWE Not(CypherRLWE C) => new(-C.a, -C.p);

    public static CypherRLWE Xor(CypherRLWE C0, CypherRLWE C1) => new(C0.a + C1.a, C0.p + C1.p);
}