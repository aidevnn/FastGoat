using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Polynoms;

public struct Fq : IGroup<EPoly<ZnInt>>
{
    public int Q { get; }
    public int P { get; }
    public int M { get; }
    public KPoly<ZnInt> F { get; }
    public EPoly<ZnInt> X { get; }

    public Fq(int q, char x = 'x')
    {
        Q = q;
        ((int p, int m), int[] coefs) = PolynomExt.GetConwayPoly(Q);
        P = p;
        M = m;
        F = new(x, ZnInt.ZnZero(p), coefs.Select(i => new ZnInt(p, i)).ToArray());
        Hash = (Q, "fq").GetHashCode();
        Name = $"F{Q}";
        FullName = Q == P ? Name : $"{Name} = F{P}({x}) = F{P}[{x}]/({F})";
        var gen = p != q ? F.X.Div(F).rem : F.One * (IntExt.Solve_k_pow_m_equal_one_mod_n_strict(p, p - 1));
        X = new(F, gen);
    }

    public IEnumerator<EPoly<ZnInt>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();

    public bool Equals(IGroup<EPoly<ZnInt>>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }
    public string FullName { get; }
    public EPoly<ZnInt> Zero => X.Zero;
    public EPoly<ZnInt> One => X.One;

    public EPoly<ZnInt> this[params ValueType[] us]
    {
        get
        {
            if (us[0] is char x && F.x.Equals(x))
                return new(F, F.X.Div(F).rem);
            if (us[0] is ValueTuple<char, int> e && F.x.Equals(e.Item1))
                return X.Pow(e.Item2);
            if (us[0] is int k)
                return X.One.Mul(k);

            throw new GroupException(GroupExceptionType.GroupDef);
        }
    }

    public IEnumerable<EPoly<ZnInt>> GetElements()
    {
        yield return Neutral();
    }

    public IEnumerable<EPoly<ZnInt>> GetGenerators()
    {
        yield return X;
    }

    public EPoly<ZnInt> Neutral() => X.One;

    public EPoly<ZnInt> Invert(EPoly<ZnInt> e) => e.Inv();

    public EPoly<ZnInt> Op(EPoly<ZnInt> e1, EPoly<ZnInt> e2) => e1 * e2;

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}