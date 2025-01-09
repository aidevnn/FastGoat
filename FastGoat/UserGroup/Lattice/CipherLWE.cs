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
        var qh = o.Mod / 2 * o;
        return new(A.Vec.Select((e, i) => i == n - 1 ? e + qh : e).ToVec());
    }

    public static CipherLWE Xor(CipherLWE A, CipherLWE B) => new(A.Vec + B.Vec);

    public static CipherLWE And(CipherLWE A, CipherLWE B, Vec<Vec<ZnInt64>> ek)
    {
        var oq = A.Vec.KOne;
        var q = oq.Mod;
        var TensAB = A.Vec.Grid2D(B.Vec).Select(e => double.Round(e.t1.Signed * e.t2.Signed * 2.0 / q)).ToArray();
        return new(ek.Select(e => (long)double.Round(e.Zip(TensAB).Sum(f => f.First.Signed * f.Second * 1.0 / q)) * oq)
            .ToVec());
    }
    
    public static CipherLWE Nand(CipherLWE A, CipherLWE B, Vec<Vec<ZnInt64>> ek) => Not(And(A, B, ek));
    public static CipherLWE Nor(CipherLWE A, CipherLWE B, Vec<Vec<ZnInt64>> ek) => And(Not(A), Not(B), ek);
    public static CipherLWE Or(CipherLWE A, CipherLWE B, Vec<Vec<ZnInt64>> ek) => Not(Nor(A, B, ek));
    public override string ToString() => $"{Vec}";
}