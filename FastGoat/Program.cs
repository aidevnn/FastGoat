using System.Numerics;
using System.Reflection;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Polynoms;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Ring.DisplayPolynomial = MonomDisplay.StarCaret;
// RecomputeAllPrimesUpTo(5000000);

void symbMinimizedForm()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var xis = Ring.Polynomial(Rational.KZero(), "x", "y", "a1", "a2", "a3", "a4", "a5", "d1", "d2", "A", "B");
    var (X, Y, a1, a2, a3, a4, a5) = Ring.EPolynomial(xis.Take(7).ToArray()).Deconstruct();
    var (d1, d2, A, B) = Ring.EPolynomial(xis.Skip(7).ToArray()).Deconstruct();
    var (x, y) = (X.Num.ExtractIndeterminate, Y.Num.ExtractIndeterminate);

    var eq1 = (X.Pow(3) + a3 * X.Pow(2) + a4 * X + a5) - (Y.Pow(2) + a1 * X * Y + a2 * Y);
    var eq2 = eq1.Substitute((Y - (a1 * X + a2) / 2, y), (X - a1.Pow(2) / 12 - a3 / 3, x));
    Console.WriteLine(eq1);
    Console.WriteLine(eq2);
    Ring.Decompose(eq2.Num, x).Item1.Println();
    // (x^3 + a3*x^2 - a1*x*y + a4*x - y^2 - a2*y + a5)
    // (1/864*a1^6 - 1/48*a1^4*x + 1/72*a1^4*a3 - 1/6*a1^2*a3*x - 1/24*a1^3*a2 + 1/18*a1^2*a3^2 + x^3 + 1/2*a1*a2*x - 1/3*a3^2*x - 1/12*a1^2*a4 - 1/6*a1*a2*a3 + 2/27*a3^3 + a4*x - y^2 + 1/4*a2^2 - 1/3*a3*a4 + a5)
    // Lines
    //     [1, 1/864*a1^6 + 1/72*a1^4*a3 - 1/24*a1^3*a2 + 1/18*a1^2*a3^2 - 1/12*a1^2*a4 - 1/6*a1*a2*a3 + 2/27*a3^3 - y^2 + 1/4*a2^2 - 1/3*a3*a4 + a5]
    //     [x, -1/48*a1^4 - 1/6*a1^2*a3 + 1/2*a1*a2 - 1/3*a3^2 + a4]
    //     [x^2, 0]
    //     [x^3, 1]
    // 

    var eq3 = (X.Pow(3) + A * X + B) - Y.Pow(2);
    var eq4 = (d1.Pow(3) * eq3.Substitute(X / d1, x)).Substitute(Y / d2, y);
    Console.WriteLine(eq3);
    Console.WriteLine(eq4);
    Console.WriteLine(eq4.Num);
    Console.WriteLine(eq4.Denom);
    Ring.Decompose(eq4.Num, x).Item1.Println();
    Console.WriteLine();
    // (-A*d1^2*d2^2*x - B*d1^3*d2^2 - d2^2*x^3 + d1^3*y^2)/(-d2^2)
    // -A*d1^2*d2^2*x - B*d1^3*d2^2 - d2^2*x^3 + d1^3*y^2
    // -d2^2
    // Lines
    //     [1, -B*d1^3*d2^2 + d1^3*y^2]
    //     [x, -A*d1^2*d2^2]
    //     [x^2, 0]
    //     [x^3, -d2^2]
    // 
}

var listElls = new List<Polynomial<Rational, Xi>>();

(EllGroup<Rational>, Func<EllPt<Rational>, EllPt<Rational>>[] revTrans)
    MinimizedForm(Polynomial<Rational, Xi> lhs, Polynomial<Rational, Xi> rhs)
{
    var F = -lhs + rhs;
    listElls.Add(F);
    var ((y, Y), (x, X)) = F.IndeterminatesAndVariables.Deconstruct();
    var ind = X.Indeterminates;
    var (xm, ym) = (new Monom<Xi>(ind, x), new Monom<Xi>(ind, y));
    var (xym, x2m) = (xm.Mul(ym), xm.Pow(2));
    var (a1, a2, a3, a4, a5) = (F[xym], F[ym], F[x2m], F[xm], F.ConstTerm);
    var A = -a1.Pow(4) / 48 - a1.Pow(2) * a3 / 6 + a1 * a2 / 2 - a3.Pow(2) / 3 + a4;
    var B = a1.Pow(6) / 864 + a1.Pow(4) * a3 / 72 - a1.Pow(3) * a2 / 24 + a1.Pow(2) * a3.Pow(2) / 18 -
        a1.Pow(2) * a4 / 12 - a1 * a2 * a3 / 6 + 2 * a3.Pow(3) / 27 + a2.Pow(2) / 4 - a3 * a4 / 3 + a5;

    var sqDivs864 = new[] { 1, 4, 9, 16, 36, 144 } // Square Divisidors of 864 
        .Select(div => (div, pow2: div * div, pow3: div * div * div)).ToArray();
    var (sqDiv, _, sqDivPow3) = sqDivs864.OrderBy(e => e.div)
        .First(div => div.pow2 % A.Denom == 0 && div.pow3 % B.Denom == 0);
    var d1 = new Rational(sqDiv);
    var d2 = new Rational(SqrtBigInt(sqDivPow3));

    Func<EllPt<Rational>, EllPt<Rational>> f1 = pt => pt.IsO ? pt : new(pt.X / d1, pt.Y / d2);
    Func<EllPt<Rational>, EllPt<Rational>> f2 = pt => pt.IsO ? pt : new(pt.X - a1.Pow(2) / 12 - a3 / 3, pt.Y);
    Func<EllPt<Rational>, EllPt<Rational>> f3 = pt => pt.IsO ? pt : new(pt.X, -pt.Y + (a1 * pt.X + a2) / 2);
    var revTrans = new[] { f1, f2, f3 };

    A *= d1.Pow(2);
    B *= d1.Pow(3);
    Console.WriteLine($"Elliptic curve      {lhs} = {rhs}");
    Console.WriteLine($"Simplified form     y^2 = {X.Pow(3) + A * X + B}");
    Console.WriteLine();

    return (new EllGroup<Rational>(A, B), revTrans);
}

void TorsionGroup((Polynomial<Rational, Xi> lhs, Polynomial<Rational, Xi> rhs) e)
{
    var F = -e.lhs + e.rhs;
    var ((y, Y), (x, X)) = F.IndeterminatesAndVariables.Deconstruct();
    var (E, revTrans) = MinimizedForm(e.lhs, e.rhs);
    var ng = EllipticCurves.NagellLutzTorsionGroup(E.A.Num, E.B.Num);
    Console.WriteLine($"TorsGroup({E}) ~ {ng.abType.Glue(" x ", "C{0}")}");

    var o = X.One;
    var torsPts = ng.gEll.ToDictionary(pt => pt, pt => revTrans.Aggregate(pt, (acc, trans) => trans(acc)));
    if (torsPts.Select(f => f.Value)
        .Any(pt => !pt.IsO && !F.Substitute(pt.X * o, x).Substitute(pt.Y * o, y).IsZero()))
        throw new();

    DisplayGroup.HeadElements(ng.gEll);
    torsPts.Println("Torsion Points");

    Console.WriteLine();
}

void testTorsionGroup()
{
    var (x, y) = Ring.Polynomial(Rational.KZero(), "x", "y").Deconstruct();

    TorsionGroup((y.Pow(2), x.Pow(3) - 2));
    TorsionGroup((y.Pow(2), x.Pow(3) + 8));
    TorsionGroup((y.Pow(2), x.Pow(3) + 4));
    TorsionGroup((y.Pow(2), x.Pow(3) + 4 * x));
    TorsionGroup((y.Pow(2) - y, x.Pow(3) - x.Pow(2)));
    TorsionGroup((y.Pow(2), x.Pow(3) + 1));
    TorsionGroup((y.Pow(2), x.Pow(3) - 43 * x + 166));
    TorsionGroup((y.Pow(2) + 7 * x * y, x.Pow(3) + 16 * x));
    TorsionGroup((y.Pow(2) + x * y + y, x.Pow(3) - x.Pow(2) - 14 * x + 29));
    TorsionGroup((y.Pow(2) + x * y, x.Pow(3) - 45 * x + 81));
    TorsionGroup((y.Pow(2) + 43 * x * y - 210 * y, x.Pow(3) - 210 * x.Pow(2)));
    TorsionGroup((y.Pow(2), x.Pow(3) - 4 * x));
    TorsionGroup((y.Pow(2) + x * y - 5 * y, x.Pow(3) - 5 * x.Pow(2)));
    TorsionGroup((y.Pow(2) + 5 * x * y - 6 * y, x.Pow(3) - 3 * x.Pow(2)));
    TorsionGroup((y.Pow(2) + 17 * x * y - 120 * y, x.Pow(3) - 60 * x.Pow(2)));
}

void testMinimizedForm()
{
    var (x, y) = Ring.Polynomial(Rational.KZero(), "x", "y").Deconstruct();

    MinimizedForm(y.Pow(2), x.Pow(3) - 2);
    MinimizedForm(y.Pow(2), x.Pow(3) + 8);
    MinimizedForm(y.Pow(2) + 5 * x * y - 6 * y, x.Pow(3) - 3 * x.Pow(2));
    MinimizedForm(y.Pow(2) + 17 * x * y - 120 * y, x.Pow(3) - 60 * x.Pow(2));
    MinimizedForm(y.Pow(2), x.Pow(3) - 4 * x - 4);
    MinimizedForm(y.Pow(2) + y, x.Pow(3) + x);
    MinimizedForm(y.Pow(2) + x * y + y, x.Pow(3) - x.Pow(2) - 29 * x - 53);
    MinimizedForm(y.Pow(2), x.Pow(3) - 11 * x - 14);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) + x.Pow(2) + x);
    MinimizedForm(y.Pow(2) + x * y + y, x.Pow(3) - 14 * x - 64);
    MinimizedForm(y.Pow(2), x.Pow(3) + 4);
    MinimizedForm(y.Pow(2) + y, x.Pow(3) + x.Pow(2) + x - 1);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) + 6 * x - 28);
    MinimizedForm(y.Pow(2), x.Pow(3) - 7 * x - 6);
    MinimizedForm(y.Pow(2), x.Pow(3) - 2 * x + 1);
    MinimizedForm(y.Pow(2) + x * y + y, x.Pow(3) + x.Pow(2));
    MinimizedForm(y.Pow(2) + x * y + y, x.Pow(3) - x.Pow(2) - 6 * x - 4);
    MinimizedForm(y.Pow(2) + y, x.Pow(3) - x.Pow(2));
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) + 15 * x + 9);
    MinimizedForm(y.Pow(2) + x * y + y, x.Pow(3) + x.Pow(2) + 1);
    MinimizedForm(y.Pow(2), x.Pow(3) + 1);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) - x.Pow(2) + 6 * x);
    MinimizedForm(y.Pow(2) + x * y + y, x.Pow(3) - 6 * x + 4);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) + 159 * x + 1737);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) - x + 137);
    MinimizedForm(y.Pow(2) + x * y + y, x.Pow(3) - x.Pow(2) - 3 * x + 3);
    MinimizedForm(y.Pow(2), x.Pow(3) + x.Pow(2) + 16 * x + 180);
    MinimizedForm(y.Pow(2), x.Pow(3) - x.Pow(2) - 4 * x + 4);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) - 34 * x + 68);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) - 4 * x - 1);
    MinimizedForm(y.Pow(2) + x * y + y, x.Pow(3) - x.Pow(2) - 14 * x + 29);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) + 108 * x + 11664);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) - 4767 * x + 127449);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) - 45 * x + 81);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) + 115 * x + 561);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) - 828 * x + 9072);
    MinimizedForm(y.Pow(2) + x * y + y, x.Pow(3) - 19 * x + 26);
    MinimizedForm(y.Pow(2) + x * y + y, x.Pow(3) - x.Pow(2) - 122 * x + 1721);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) - 361 * x + 2585);
    MinimizedForm(y.Pow(2) + x * y + y, x.Pow(3) + 1922 * x + 20756);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) - 1070 * x + 7812);
    MinimizedForm(y.Pow(2) + x * y, x.Pow(3) - 8696090 * x + "9838496100");
}

{
    symbMinimizedForm();
    testMinimizedForm();
    testTorsionGroup();

    Console.WriteLine();

    // export to SageMath
    foreach (var (idx, F) in listElls.Index())
    {
        Console.WriteLine($"E{idx:00}=EllipticCurve({F})");
        Console.WriteLine($"E{idx:00}.torsion_subgroup()");
        Console.WriteLine($"E{idx:00}.short_weierstrass_model()");
        Console.WriteLine($"E{idx:00}.torsion_points()");
        Console.WriteLine();
    }
}