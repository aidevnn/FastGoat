using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.UserGroup.Lattice;

public struct CipherLWE
{
    public Vec<ZnInt64> Vec { get; }

    public CipherLWE(Vec<ZnInt64> vec)
    {
        Vec = vec;
    }

    public static CipherLWE Not(CipherLWE A)
    {
        var n = A.Vec.Length;
        var o = A.Vec.KOne;
        var qh = o.P / 2 * o;
        return new(A.Vec.Select((e, i) => i == n - 1 ? e + qh : e).ToVec());
    }

    public static CipherLWE Xor(CipherLWE A, CipherLWE B) => new(A.Vec + B.Vec);

    public static CipherLWE And(CipherLWE A, CipherLWE B, Vec<Vec<ZnInt64>> ek)
    {
        var oq = A.Vec.KOne;
        var oq2 = ek[0].KOne;
        var q = oq.P;
        var cstar = A.Vec.Grid2D(B.Vec).Select(e => LWE.Signed(e.t1) * LWE.Signed(e.t2) * oq2).ToVec();
        return new((cstar * ek).Select(e => (LWE.Signed(e.Sum()) * 2 / q) * oq).ToVec());
    }

    public static CipherLWE Nand(CipherLWE A, CipherLWE B, Vec<Vec<ZnInt64>> ek) => Not(And(A, B, ek));
    
    public override string ToString() => $"{Vec}";
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