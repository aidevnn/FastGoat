using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Polynoms;

public class GFp : IGroup<FpPolynom>
{
    public GFp(string name, (char x, int p) xp, int[] polynom)
    {
        if (!IntExt.Primes10000.Contains(xp.p))
            throw new GroupException(GroupExceptionType.GroupDef);

        P = xp.p;
        X = xp.x;
        Name = name;

        polynom = PolynomExt.TrimPoly(polynom);
        if (polynom.Length < 2 || polynom.Last() != 1 ||
            polynom.Any(e => IntExt.AmodP(e, P) != e))
            throw new GroupException(GroupExceptionType.NotIrreductiblePolynom);

        var seq = Enumerable.Range(1, P - 1).Select(a => polynom.Select((e, i) => e * (int)Math.Pow(a, i)).Sum())
            .ToArray();

        Poly = new FpPolynom(P, X, polynom);
        M = polynom.Length - 1;

        Hash = (Poly.Hash, xp.x).GetHashCode();
        _neutral = FpPolynom.One(P, X);

        Elements = Group.GenerateElements(this, FpPolynom.Fpx(P, X)).ToHashSet();
    }

    public FpPolynom Poly { get; }

    private char X { get; }
    private readonly FpPolynom _neutral;
    private HashSet<FpPolynom> Elements { get; }
    public int M { get; }
    public int P { get; }

    public IEnumerator<FpPolynom> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<FpPolynom>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public FpPolynom this[params ValueType[] us]
    {
        get
        {
            var intUs = us.Select(Convert.ToInt32).ToArray();
            return new FpPolynom(P, X, intUs).DivideBy(Poly).rem;
        }
    }

    public IEnumerable<FpPolynom> GetElements() => Elements;

    public IEnumerable<FpPolynom> GetGenerators()
    {
        yield return FpPolynom.Fpx(P, X);
    }

    public FpPolynom Neutral() => _neutral;

    public FpPolynom Invert(FpPolynom e)
    {
        var (x, y) = FpPolynom.Bezout(e, Poly);
        var gcd = e * x + y * Poly;
        return (x / gcd).DivideBy(Poly).rem;
    }

    public FpPolynom Op(FpPolynom e1, FpPolynom e2)
    {
        return (e1 * e2).DivideBy(Poly).rem;
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}