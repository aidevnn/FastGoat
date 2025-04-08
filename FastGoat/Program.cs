using System.Collections;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void symbWeierstrassForm()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var xis = Ring.Polynomial(Rational.KZero(), MonomOrder.GrLex,
        "a1", "a2", "a3", "a4", "a5", "d1", "d2", "A", "B", "C", "x", "y");
    var (a1, a2, a3, a4, a5, d1, d2) = Ring.EPolynomial(xis.Take(7).ToArray()).Deconstruct();
    var (A, B, C, X, Y) = Ring.EPolynomial(xis.Skip(7).ToArray()).Deconstruct();
    var (x, y) = (X.Num.ExtractIndeterminate, Y.Num.ExtractIndeterminate);

    var eqEll = (X.Pow(3) + a2 * X.Pow(2) + a4 * X + a5) - (Y.Pow(2) + a1 * X * Y + a3 * Y);
    
    var eqLong1 = eqEll.Substitute(Y - (a1 * X + a3) / 2, y);
    Console.WriteLine(eqEll);
    Console.WriteLine(eqLong1);
    Ring.Decompose(eqLong1.Num, x).Item1.Println();
    // (-a1*x*y + a2*x^2 + x^3 - a3*y + a4*x - y^2 + a5)
    // (1/4*a1^2*x^2 + 1/2*a1*a3*x + a2*x^2 + x^3 - 1/4*a2^2 + 1/2*a2*a3 + a2*y - a3*y + a4*x - y^2 + a5)
    // Lines
    //     [1, -1/4*a2^2 + 1/2*a2*a3 + a2*y - a3*y - y^2 + a5]
    //     [x, 1/2*a1*a3 + a4]
    //     [x^2, 1/4*a1^2 + a2]
    //     [x^3, 1]
    // 

    var eqLong2 = (X.Pow(3) + A * X * X + B * X + C) - Y.Pow(2);
    var eqLong3 = (d1.Pow(3) * eqLong2.Substitute(X / d1, x)).Substitute(Y / d2, y);
    Console.WriteLine(eqLong2);
    Console.WriteLine(eqLong3);
    Console.WriteLine(eqLong3.Num);
    Console.WriteLine(eqLong3.Denom);
    Ring.Decompose(eqLong3.Num, x).Item1.Println();
    Console.WriteLine();
    // (A*x^2 + x^3 + B*x - y^2 + C)
    // (-C*d1^3*d2^2 - B*d1^2*d2^2*x - A*d1*d2^2*x^2 + d1^3*y^2 - d2^2*x^3)/(-d2^2)
    // -C*d1^3*d2^2 - B*d1^2*d2^2*x - A*d1*d2^2*x^2 + d1^3*y^2 - d2^2*x^3
    // -d2^2
    // Lines
    //     [1, -C*d1^3*d2^2 + d1^3*y^2]
    //     [x, -B*d1^2*d2^2]
    //     [x^2, -A*d1*d2^2]
    //     [x^3, -d2^2]
    // 
    
    var eqShort1 = eqEll.Substitute((Y - (a1 * X + a3) / 2, y), (X - a1.Pow(2) / 12 - a2 / 3, x));
    Console.WriteLine(eqEll);
    Console.WriteLine(eqShort1);
    Ring.Decompose(eqShort1.Num, x).Item1.Println();
    // (-a1*x*y + a2*x^2 + x^3 - a3*y + a4*x - y^2 + a5)
    // (1/864*a1^6 + 1/72*a1^4*a2 - 1/48*a1^4*x - 1/24*a1^3*a3 + 1/18*a1^2*a2^2 - 1/6*a1^2*a2*x - 1/12*a1^2*a4 - 1/6*a1*a2*a3 + 1/2*a1*a3*x + 2/27*a2^3 - 1/3*a2^2*x + x^3 - 1/3*a2*a4 + 1/4*a3^2 + a4*x - y^2 + a5)
    // Lines
    //     [1, 1/864*a1^6 + 1/72*a1^4*a2 - 1/24*a1^3*a3 + 1/18*a1^2*a2^2 - 1/12*a1^2*a4 - 1/6*a1*a2*a3 + 2/27*a2^3 - 1/3*a2*a4 + 1/4*a3^2 - y^2 + a5]
    //     [x, -1/48*a1^4 - 1/6*a1^2*a2 + 1/2*a1*a3 - 1/3*a2^2 + a4]
    //     [x^2, 0]
    //     [x^3, 1]
    // 

    var eqShort2 = (X.Pow(3) + A * X + B) - Y.Pow(2);
    var eqShort3 = (d1.Pow(3) * eqShort2.Substitute(X / d1, x)).Substitute(Y / d2, y);
    Console.WriteLine(eqShort2);
    Console.WriteLine(eqShort3);
    Console.WriteLine(eqShort3.Num);
    Console.WriteLine(eqShort3.Denom);
    Ring.Decompose(eqShort3.Num, x).Item1.Println();
    // (-B*d1^3*d2^2 - A*d1^2*d2^2*x + d1^3*y^2 - d2^2*x^3)/(-d2^2)
    // -B*d1^3*d2^2 - A*d1^2*d2^2*x + d1^3*y^2 - d2^2*x^3
    // -d2^2
    // Lines
    //     [1, -B*d1^3*d2^2 + d1^3*y^2]
    //     [x, -A*d1^2*d2^2]
    //     [x^2, 0]
    //     [x^3, -d2^2]
    // 

    Console.WriteLine();
}

ConcreteGroup<EllPt<ZnInt>> EllFp(EllGroupLong<ZnInt> E)
{
    var (A, B) = (E.ShortForm.A, E.ShortForm.B);
    var p = A.P;
    var ell = p.Range().Select(k => new ZnInt(p, k))
        .Select(x => (x, y2: x.Pow(3) + A * x + B))
        .Select(e => (e.x, y: NumberTheory.SqrtModANTV1(e.y2.K, p) * e.x.One))
        .Select(e => E.ConvertFromShort(new(e.x, e.y)))
        .Where(e => E.Contains(e.X, e.Y))
        .Order()
        .ToArray();

    return Group.Generate(E, ell);
}

ConcreteGroup<EllPt<ZnInt>> EllFpLongForm(int p, int[] curve)
{
    var (a1, a2, a3, a4, a5) = curve.Select(i => new ZnInt(p, i)).Deconstruct();
    var E = new EllGroupLong<ZnInt>(a1, a2, a3, a4, a5);

    var gEll = EllFp(E);
    var abType = Group.AbelianGroupType(gEll);
    Console.WriteLine($"{gEll} ~ {abType.Glue(" x ", "C{0}")}");
    DisplayGroup.Generators(gEll, showBaseGroup: false);
    return gEll;
}

IEnumerable<EllPt<Rational>> SolveIntegralPoints(EllGroupLong<Rational> E, bool show = false)
{
    var disc = E.Disc;
    var (A, B, C, _, _) = E.LongForm;
    Console.WriteLine($"{E} ~ {E.LongFormStr} ~ {E.ShortFormStr}");
    // var (a1, a2, a3, a4, a5) = E.Coefs;
    // var (X, Y) = Ring.Polynomial(Rational.KOne(), "X", "Y").Deconstruct();
    // var F = X.Pow(3) + a2 * X.Pow(2) + a4 * X + a5 - (Y.Pow(2) + a1 * X * Y + a3 * Y);
    // var disc1 = Ring.Discriminant(Ring.Discriminant(F, Y) / 16, X).ConstTerm / 16;
    // Console.WriteLine(new { disc, disc1, div = disc / disc1 });
    var r = IntExt.PrimesDec(BigInteger.Abs(disc.Num))
        .Aggregate(BigInteger.One, (acc, r) => acc * BigInteger.Pow(r.Key, r.Value / 2 + r.Value % 2));
    var divs = IntExt.DividorsBigInt(16 * r).Where(y => (256 * disc.Num) % (y * y) == 0).Order()
        .Select(y => new Rational(y)).ToArray();

    var x = FG.QPoly();
    foreach (var y in divs.Prepend("0"))
    {
        var P = x.Pow(3) + A * x * x + B * x + C;
        var sols = IntFactorisation.FactorsQ(P - y.Pow(2));
        if (show)
            sols.Println($"Y = {y}, solve {y.Pow(2)} = {P}");

        var ellpts = sols.Where(e => e.Item1.Degree == 1).Select(e => new EllPt<Rational>(-e.Item1[0], y));
        foreach (var pt in ellpts)
        {
            yield return E.ConvertFromLong(pt);
            yield return E.ConvertFromLong(new(pt.X, -pt.Y));
        }
    }
}

(ConcreteGroup<EllPt<Rational>> gEll, HashSet<EllPt<Rational>> pts) NagellLutzTorsionGroup(BigInteger[] curve)
{
    var (a1, a2, a3, a4, a5) = curve.Select(e => new Rational(e)).Deconstruct();
    var E = new EllGroupLong<Rational>(a1, a2, a3, a4, a5);
    var ellpts = SolveIntegralPoints(E).ToArray();
    var pts = ellpts.ToHashSet();
    var set = new List<EllPt<Rational>>() { E.O };

    foreach (var pt in ellpts.Where(pt => E.ConvertToShort(pt).IsIntegral()))
    {
        var acc = pt;
        for (int i = 1; i <= 4; i++)
        {
            acc = E.Times(acc, 2);
            if (acc.IsO || !E.ConvertToShort(acc).IsIntegral())
                break;
        }

        if (E.ConvertToShort(acc).IsIntegral())
            set.Add(pt);
    }

    var gEll = Group.Generate(E, set.ToArray());
    var abType = Group.AbelianGroupType(gEll);
    DisplayGroup.HeadElements(gEll);
    Console.WriteLine($"Tors({gEll}) ~ {abType.Glue(" x ", "C{0}")}");
    Console.WriteLine();

    return (gEll, pts);
}

void testEllFp(int[] curve)
{
    var (a1, a2, a3, a4, a5) = curve.Select(i => new Rational(i)).Deconstruct();
    var E0 = new EllGroupLong<Rational>(a1, a2, a3, a4, a5);
    foreach (var p in Primes10000.Where(p => p > 3 && E0.Disc.Num % p != 0).Take(10))
    {
        var gEll = EllFpLongForm(p, curve);
        var E = (EllGroupLong<ZnInt>)gEll.BaseGroup;

        var gEll2 = EllFpLongForm(p, [0, 0, 0, E.ShortForm.A.K, E.ShortForm.B.K]);
        Console.WriteLine();

        if (gEll.Any(pt => gEll.ElementsOrders[pt] != gEll2.ElementsOrders[E.ConvertToShort(pt)]))
            throw new();
    }
}

void testEllTorsNG()
{
    GlobalStopWatch.AddLap();

    // Torsion C1
    NagellLutzTorsionGroup([0, 0, 0, -4, -4]);
    NagellLutzTorsionGroup([0, 0, 1, 1, 0]);
    NagellLutzTorsionGroup([1, -1, 1, -29, -53]);

    // Torsion C2
    NagellLutzTorsionGroup([0, 0, 0, -11, -14]);
    NagellLutzTorsionGroup([1, 1, 0, 1, 0]);
    NagellLutzTorsionGroup([1, 0, 1, -14, -64]);

    // Torsion C3
    NagellLutzTorsionGroup([0, 0, 0, 0, 4]);
    NagellLutzTorsionGroup([0, 1, 1, 1, -1]);
    NagellLutzTorsionGroup([1, 0, 0, 6, -28]);

    // Torsion C4, C2 x C2
    NagellLutzTorsionGroup([0, 0, 0, -7, -6]);
    NagellLutzTorsionGroup([0, 0, 0, -2, 1]);
    NagellLutzTorsionGroup([1, 1, 1, 0, 0]);
    NagellLutzTorsionGroup([1, -1, 1, -6, -4]);

    // Torsion C5
    NagellLutzTorsionGroup([0, -1, 1, 0, 0]);
    NagellLutzTorsionGroup([1, 0, 0, 15, 9]);
    NagellLutzTorsionGroup([1, 1, 1, 0, 1]);

    // Torsion C6
    NagellLutzTorsionGroup([0, 0, 0, 0, 1]);
    NagellLutzTorsionGroup([1, -1, 0, 6, 0]);
    NagellLutzTorsionGroup([1, 0, 1, -6, 4]);

    // Torsion C7
    NagellLutzTorsionGroup([1, 0, 0, 159, 1737]);
    NagellLutzTorsionGroup([1, 0, 0, -1, 137]);
    NagellLutzTorsionGroup([1, -1, 1, -3, 3]);

    // Torsion C8, C4 x C2
    NagellLutzTorsionGroup([0, 1, 0, 16, 180]);
    NagellLutzTorsionGroup([0, -1, 0, -4, 4]);
    NagellLutzTorsionGroup([1, 0, 0, -34, 68]);
    NagellLutzTorsionGroup([1, 0, 0, -4, -1]);

    // Torsion C9
    NagellLutzTorsionGroup([1, -1, 1, -14, 29]);
    NagellLutzTorsionGroup([1, 0, 0, 108, 11664]);
    NagellLutzTorsionGroup([1, 0, 0, -4767, 127449]);

    // Torsion C10
    NagellLutzTorsionGroup([1, 0, 0, -45, 81]);
    NagellLutzTorsionGroup([1, 0, 0, 115, 561]);
    NagellLutzTorsionGroup([1, 0, 0, -828, 9072]);

    // Torsion C12, C6 x C2
    NagellLutzTorsionGroup([1, 0, 1, -19, 26]);
    NagellLutzTorsionGroup([1, -1, 1, -122, 1721]);
    NagellLutzTorsionGroup([1, 0, 0, -361, 2585]);
    NagellLutzTorsionGroup([1, 0, 1, 1922, 20756]);

    // Torsion C8 x C2
    NagellLutzTorsionGroup([1, 0, 0, -1070, 7812]);
    NagellLutzTorsionGroup([1, 0, 0, -8696090, 9838496100]);

    GlobalStopWatch.Show(); // #  Time:5.955s
    Console.WriteLine();
}

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    GlobalStopWatch.Restart();
    symbWeierstrassForm();

    EllFpLongForm(p: 7, curve: [1, 1, 0, 1, 1]);
    EllFpLongForm(p: 17, curve: [1, 1, 0, 1, 1]);
    testEllFp(curve: [0, 1, 0, 1, 3]);
    testEllFp(curve: [1, -1, 1, -29, -53]);
    
    NagellLutzTorsionGroup([0, 0, 0, -36, 0]);
    NagellLutzTorsionGroup([0, 13, 0, -240, 576]);
    NagellLutzTorsionGroup([5, -3, -6, 0, 0]);

    testEllTorsNG();
    testEllTorsNG();
    testEllTorsNG();
}

public struct EllGroupLong<T> : IGroup<EllPt<T>> where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>
{
    public EllGroupLong(T a1, T a2, T a3, T a4, T a5)
    {
        var A1 = -a1.Pow(4) / 48 - a1.Pow(2) * a2 / 6 + a1 * a3 / 2 - a2.Pow(2) / 3 + a4;
        var B1 = a1.Pow(6) / 864 + a1.Pow(4) * a2 / 72 - a1.Pow(3) * a3 / 24 + a1.Pow(2) * a2.Pow(2) / 18 -
            a1.Pow(2) * a4 / 12 - a1 * a3 * a2 / 6 + 2 * a2.Pow(3) / 27 + a3.Pow(2) / 4 - a2 * a4 / 3 + a5;

        Disc = -4 * A1.Pow(3) - 27 * B1.Pow(2);
        if (Disc.IsZero())
            throw new GroupException(GroupExceptionType.GroupDef);

        T d11 = a1.One, d12 = a1.One;
        if (A1 is Rational _A1 && B1 is Rational _B1)
        {
            var sqDivs864 = new[] { 1, 4, 9, 16, 36, 144 } // Square Divisidors of 864 
                .Select(div => (div, pow2: div * div, pow3: div * div * div)).ToArray();
            var (sqDiv, _, sqDivPow3) = sqDivs864.OrderBy(f => f.div)
                .First(div => div.pow2 % _A1.Denom == 0 && div.pow3 % _B1.Denom == 0);
            dynamic _d11 = new Rational(sqDiv);
            dynamic _d12 = new Rational(IntExt.SqrtBigInt(sqDivPow3));
            (d11, d12) = (_d11, _d12);
        }

        A1 *= d11.Pow(2);
        B1 *= d11.Pow(3);

        var (A2, B2, C2) = (a1.Pow(2) / 4 + a2, a1 * a3 / 2 + a4, a3.Pow(2) / 4 + a5);

        T d21 = a1.One, d22 = a1.One;
        if (A2 is Rational A0 && B2 is Rational B0 && C2 is Rational C0 &&
            (!A0.IsInteger() || !B0.IsInteger() || !C0.IsInteger()))
            (d21, d22) = (4 * a1.One, 8 * a1.One);

        A2 *= d21;
        B2 *= d21.Pow(2);
        C2 *= d21.Pow(3);

        Coefs = (a1, a2, a3, a4, a5);
        ShortForm = (A1, B1, d11, d12);
        LongForm = (A2, B2, C2, d21, d22);

        Field = typeof(T).Name;
        if (A1 is Rational)
            Field = "Q";
        else if (A1 is ZnInt _a1)
            Field = $"Z/{_a1.P}Z";
        else if (A1 is ZnBigInt _a2)
            Field = $"Z/{_a2.Mod}Z";

        Name = $"Ell[{a1},{a2},{a3},{a4},{a5}]({Field})".Replace(" ", "");
        Hash = (Coefs, ShortForm).GetHashCode();
    }

    public EllGroupLong(T a, T b) : this(a.Zero, a.Zero, a.Zero, a, b)
    {
    }

    public EllGroupLong(T a, T b, T c) : this(a.Zero, a, a.Zero, b, c)
    {
    }

    public IEnumerator<EllPt<T>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<EllPt<T>>? other) => other?.Hash == Hash;

    public EllPt<T> ConvertToShort(EllPt<T> pt)
    {
        if (pt.IsO)
            return pt;

        var _x = (pt.X + a1.Pow(2) / 12 + a2 / 3) * ShortForm.d1;
        var _y = (pt.Y + (a1 * pt.X + a3) / 2) * ShortForm.d2;
        return new(_x, _y);
    }

    public EllPt<T> ConvertFromShort(EllPt<T> pt)
    {
        if (pt.IsO)
            return pt;

        var _x = pt.X / ShortForm.d1 - a1.Pow(2) / 12 - a2 / 3;
        var _y = pt.Y / ShortForm.d2;
        return new(_x, _y - (a1 * _x + a3) / 2);
    }

    public EllPt<T> ConvertToLong(EllPt<T> pt)
    {
        if (pt.IsO)
            return pt;

        return new(pt.X * LongForm.d1, (pt.Y + (a1 * pt.X + a3) / 2) * LongForm.d2);
    }

    public EllPt<T> ConvertFromLong(EllPt<T> pt)
    {
        if (pt.IsO)
            return pt;

        return new(pt.X / LongForm.d1, (pt.Y / LongForm.d2 - (a1 * pt.X / LongForm.d1 + a3) / 2));
    }

    public int Hash { get; }
    public string Name { get; }
    public string Field { get; }
    public T Disc { get; }
    public (T a1, T a2, T a3, T a4, T a5) Coefs { get; }
    public (T A, T B, T d1, T d2) ShortForm { get; }
    public (T A, T B, T C, T d1, T d2) LongForm { get; }
    public string ShortFormStr => $"Ell[{ShortForm.A},{ShortForm.B}]({Field})".Replace(" ", "");
    public string LongFormStr => $"Ell[0,{LongForm.A},0,{LongForm.B},{LongForm.C}]({Field})".Replace(" ", "");
    public T a1 => Coefs.a1;
    public T a2 => Coefs.a2;
    public T a3 => Coefs.a3;
    public T a4 => Coefs.a4;
    public T a5 => Coefs.a5;

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

    public bool Contains(T X, T Y)
    {
        var lhs = Y * Y + a1 * X * Y + a3 * Y;
        var rhs = X.Pow(3) + a2 * X * X + a4 * X + a5;
        return lhs.Equals(rhs);
    }

    public EllPt<T> O => new();
    public EllPt<T> Neutral() => O;

    public EllPt<T> Invert(EllPt<T> P)
    {
        if (P.IsO)
            return P;

        if (!Contains(P.X, P.Y))
            throw new GroupException(GroupExceptionType.GroupDef);

        var P1 = ConvertToShort(P);
        return ConvertFromShort(new EllPt<T>(P1.X, -P1.Y));
    }

    public EllPt<T> Op(EllPt<T> e1, EllPt<T> e2)
    {
        if (e1.IsO)
            return e2;

        if (e2.IsO)
            return e1;

        if (!Contains(e1.X, e1.Y) || !Contains(e2.X, e2.Y))
        {
            Console.WriteLine(new { e1, e2, E = this });
            throw new GroupException(GroupExceptionType.GroupDef);
        }

        var (e1_, e2_) = (ConvertToShort(e1), ConvertToShort(e2));
        var (x1, y1, x2, y2) = (e1_.X, e1_.Y, e2_.X, e2_.Y);
        if (!x1.Equals(x2))
        {
            var alpha = (y2 - y1) / (x2 - x1);
            var x3 = alpha.Pow(2) - x1 - x2;
            var y3 = -y1 + alpha * (x1 - x3);
            return ConvertFromShort(new(x3, y3));
        }
        else
        {
            if (!y1.Equals(y2.Opp()))
            {
                var alpha = (3 * x1.Pow(2) + ShortForm.A) / (2 * y1);
                var x3 = alpha.Pow(2) - 2 * x1;
                var y3 = -y1 + alpha * (x1 - x3);
                return ConvertFromShort(new(x3, y3));
            }
            else
                return new();
        }
    }

    public override int GetHashCode() => Hash;

    public override string ToString() => Name;
}