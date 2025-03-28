using System.Collections;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Padic;

namespace FastGoat.UserGroup.EllCurve;

// Daniel Guin - Thomas Hausberger, Algebre Tome 1, page 166
public struct EllGroup<T> : IGroup<EllPt<T>> where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>
{
    public EllGroup(T a, T b)
    {
        var disc = 4 * a.Pow(3) + 27 * b.Pow(2);
        if (disc.IsZero())
            throw new GroupException(GroupExceptionType.GroupDef);
        
        (A, B) = (a, b);
        var field = typeof(T).Name;
        if (a is Rational a0)
            field = "Q";
        else if (a is ZnInt a1)
            field = $"Z/{a1.P}Z";
        else if (a is ZnBigInt a2)
            field = $"Z/{a2.Mod}Z";
        
        Name = $"Ell[{a},{b}]({field})";
        Hash = (a, b).GetHashCode();
    }
    public IEnumerator<EllPt<T>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<EllPt<T>>? other) => other?.Hash == Hash;
    
    public int Hash { get; }
    public string Name { get; }
    public T A { get; }
    public T B { get; }

    public EllPt<T> this[params ValueType[] us]
    {
        get
        {
            if (us.Length == 2 && us[0] is T x && us[1] is T y)
            {
                if (!Contains(x, y))
                    throw new GroupException(GroupExceptionType.GroupDef);
                
                return new EllPt<T>(x, y);
            }
            
            throw new GroupException(GroupExceptionType.GroupDef);
        }
    }

    public IEnumerable<EllPt<T>> GetElements()
    {
        yield return new EllPt<T>();
    }

    public IEnumerable<EllPt<T>> GetGenerators()
    {
        yield return new EllPt<T>();
    }

    public bool Contains(T X, T Y) => (Y * Y).Equals(X.Pow(3) + A * X + B);
    public EllPt<T> O => new();
    public EllPt<T> Neutral() => O;

    public EllPt<T> Invert(EllPt<T> P)
    {
        if (P.IsO)
            return P;
        
        if (!Contains(P.X, P.Y))
            throw new GroupException(GroupExceptionType.GroupDef);
        
        return new EllPt<T>(P.X, -P.Y);
    }

    public EllPt<T> Op(EllPt<T> e1, EllPt<T> e2)
    {
        if (e1.IsO)
            return e2;

        if (e2.IsO)
            return e1;

        var (x1, y1, x2, y2) = (e1.X, e1.Y, e2.X, e2.Y);
        if (!Contains(x1, y1) || !Contains(x2, y2))
        {
            Console.WriteLine(new { e1, e2, E = this });
            throw new GroupException(GroupExceptionType.GroupDef);
        }
        
        if (!x1.Equals(x2))
        {
            var alpha = (y2 - y1) / (x2 - x1);
            var x3 = alpha.Pow(2) - x1 - x2;
            var y3 = -y1 + alpha * (x1 - x3);
            return new(x3, y3);
        }
        else
        {
            if (!y1.Equals(y2.Opp()))
            {
                var alpha = (3 * x1.Pow(2) + A) / (2 * y1);
                var x3 = alpha.Pow(2) - 2 * x1;
                var y3 = -y1 + alpha * (x1 - x3);
                return new(x3, y3);
            }
            else
                return new();
        }
    }

    public override int GetHashCode() => Hash;

    public override string ToString() => Name;
}