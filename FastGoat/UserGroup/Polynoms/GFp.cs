using System.Collections;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Padic;

namespace FastGoat.UserGroup.Polynoms;

public class GFp : IGroup<EPoly<ZnInt>>
{
    public GFp(string name, EPoly<ZnInt> e)
    {
        if (!IntExt.Primes10000.Contains(e.P))
            throw new GroupException(GroupExceptionType.GroupDef);

        P = e.P;
        Name = name;

        F = e.F;
        X = e.X;
        M = F.Degree;

        Hash = e.Hash;
    }

    public GFp(string name, KPoly<ZnInt> f) : this(name, FG.EPoly(f))
    {
    }

    public GFp(KPoly<ZnInt> f) : this($"F{f.P}[{f.x}]/({f})", FG.EPoly(f))
    {
    }

    public GFp(EPoly<ZnInt> e) : this($"F{e.P}[{e.F.x}]/({e.F})", e)
    {
    }

    public KPoly<ZnInt> F { get; }
    public EPoly<ZnInt> X { get; }
    public int M { get; }
    public int P { get; }

    public IEnumerator<EPoly<ZnInt>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<EPoly<ZnInt>>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public EPoly<ZnInt> this[params ValueType[] us]
    {
        get
        {
            if (us[0] is char x && F.x.Equals(x))
                return X;
            if (us[0] is EPoly<ZnInt> e && e.P == P && e.Poly.Equals(F) && !e.IsZero())
                return e;
            if (us[0] is int k && k != 0)
                return X * k;

            throw new GroupException(GroupExceptionType.GroupDef);
        }
    }

    public IEnumerable<EPoly<ZnInt>> GetElements()
    {
        yield return Neutral();
    }

    public IEnumerable<EPoly<ZnInt>> GetGenerators()
    {
        yield return NumberTheory.PrimitiveRoot(X);
    }

    public EPoly<ZnInt> Neutral() => X.One;

    public EPoly<ZnInt> Invert(EPoly<ZnInt> e) => e.Inv();

    public EPoly<ZnInt> Op(EPoly<ZnInt> e1, EPoly<ZnInt> e2) => e1 * e2;

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}