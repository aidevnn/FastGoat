using System.Numerics;
using System.Text.Json.Serialization;
using System.Text.RegularExpressions;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;
using Group = FastGoat.Structures.Group;
using RegXGroup = System.Text.RegularExpressions.Group;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

(KPoly<Rational> X2, KPoly<Rational> Y2) EllP2(EllCoefs<Rational> E)
{
    var (a1, a2, a3, a4, a6) = E.Model;
    var xis = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "Y", "X");
    var (Y0, X0) = xis.Deconstruct();
    var basis = new PolynomialBasis<Rational, Xi>((Y0.Pow(2) + a1 * X0 * Y0 + a3 * Y0) -
                                                  (X0.Pow(3) + a2 * X0.Pow(2) + a4 * X0 + a6));
    
    var (Y, X) = Ring.EPolynomial(xis).Deconstruct();
    var Ell = new EllGroup<EPolynomial<Rational>>(a1 * X.One, a2 * X.One, a3 * X.One, a4 * X.One, a6 * X.One);
    var P = new EllPt<EPolynomial<Rational>>(X, Y);

    var P1 = Ell.ConvertToShort(P);
    var (x1, y1) = P1;
    var alpha = (3 * x1.Pow(2) + Ell.ShortForm.A) / (2 * y1);
    var x2 = alpha.Pow(2) - 2 * x1;
    var y2 = -y1 + alpha * (x1 - x2);
    var P2 = Ell.ConvertFromShort(new(x2, y2));
    var xn = basis.Rem(P2.X.Num);
    var yn = basis.Rem(P2.Y.Num);
    P2 = new(new(xn, P2.X.Denom), new(yn.Monic(), P2.Y.Denom / yn.LeadingDetails.lc));
    Console.WriteLine($"{E.Eq}, P = {P} and 2P = {P2}");
    return (xn.ToKPoly("X").Monic, yn.ToKPoly("Y").Monic);
}

KPoly<KPoly<Rational>> Simplify(KPoly<KPoly<Rational>> P, KPoly<KPoly<Rational>> R)
{
    if (P.Degree < 2)
        return P;

    var (quo, rem) = P.Div(P.X.Pow(2));
    return Simplify(quo * R + rem, R);
}

(KPoly<KPoly<Rational>> R, Dictionary<int, KPoly<KPoly<Rational>>> psi, Dictionary<int, KPoly<Rational>> f) 
    DivisionPolynomial(EllGroup<Rational> E, int nmax)
{
    var x = FG.QPoly('X');
    var Y = FG.KPoly('Y', x);
    var X = x * Y.One;
    var (A, B) = (E.ShortForm.A * x.One * Y.One, E.ShortForm.B * x.One * Y.One);
    var R = X.Pow(3) + A * X + B;

    var psi = new Dictionary<int, KPoly<KPoly<Rational>>>();
    (psi[0], psi[1], psi[2]) = (Y.Zero, Y.One, 2 * Y);
    psi[3] = 3 * X.Pow(4) + 6 * A * X.Pow(2) + 12 * B * X - A.Pow(2);
    psi[4] = 4 * Y * (X.Pow(6) + 5 * A * X.Pow(4) + 20 * B * X.Pow(3) - 5 * A.Pow(2) * X.Pow(2) - 4 * A * B * X -
                      8 * B.Pow(2) - A.Pow(3));

    for (int n = 2; n < nmax / 2; n++)
    {
        if (n > 2)
            psi[2 * n] = Simplify(psi[n] * (psi[n + 2] * psi[n - 1].Pow(2) - psi[n - 2] * psi[n + 1].Pow(2)) / (2 * Y),
                R);
        psi[2 * n + 1] = Simplify(psi[n + 2] * psi[n].Pow(3) - psi[n + 1].Pow(3) * psi[n - 1], R);
    }

    var f = psi.ToDictionary(e => e.Key, e => e.Key % 2 == 0 ? (e.Value / Y)[0].Clone() : e.Value[0].Clone());
    return (R, psi, f);
}

void nP(int n, Dictionary<int, KPoly<KPoly<Rational>>> psi, KPoly<KPoly<Rational>> R)
{
    var X = FG.KFracPoly(R.KOne.X).X;
    var Y = FG.KFracPoly('Y', X).X;

    var tmp1 = Simplify(psi[n - 1] * psi[n + 1], R).ToFrac(Y) / Simplify(psi[n].Pow(2), R).ToFrac(Y);
    var nPX = X - tmp1;

    var denom = n % 2 == 0
        ? Y * Simplify(4 * psi[n].Pow(3), R).ToFrac(Y)
        : Simplify(4 * R.X * psi[n].Pow(3), R).ToFrac(Y);
    var tmp2 = Simplify(psi[n + 2] * psi[n - 1].Pow(2) - psi[n - 2] * psi[n + 1].Pow(2), R).ToFrac(Y);
    var nPY = tmp2 / denom;

    Console.WriteLine($"P=({X},{Y}) and {n}P=({nPX},{nPY})");
    Console.WriteLine();
}

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var E = EC.EllCoefs([1, 0]);
    Console.WriteLine($"Ell({E.ModelStr})(Q) {E.Eq}");
    var (R, psi, divPolys) = DivisionPolynomial(E.ToEllGroup(), 8);
    psi.Println("psi");
    divPolys.Println("divPolys");
    divPolys.Where(e => e.Key != 0)
        .ToDictionary(e => e.Key, e => (e.Value.Degree, (e.Key.Pow(2) - 1 - 3 * ((e.Key + 1) % 2)) / 2))
        .Println("Degree of divPolys");

    EllP2(E);
    nP(2, psi, R);
    nP(3, psi, R);
}

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var E = EC.EllCoefs([1, 1]);
    Console.WriteLine($"Ell({E.ModelStr})(Q) {E.Eq}");
    var (R, psi, divPolys) = DivisionPolynomial(E.ToEllGroup(), 8);
    psi.Println("psi");
    divPolys.Println("divPolys");
    divPolys.Where(e => e.Key != 0)
        .ToDictionary(e => e.Key, e => (e.Value.Degree, (e.Key.Pow(2) - 1 - 3 * ((e.Key + 1) % 2)) / 2))
        .Println("Degree of divPolys");
    
    EllP2(E);
    nP(2, psi, R);
    nP(3, psi, R);
}