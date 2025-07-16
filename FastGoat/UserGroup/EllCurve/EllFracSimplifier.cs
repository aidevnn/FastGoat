using FastGoat.Structures;
using FastGoat.Structures.VecSpace;

namespace FastGoat.UserGroup.EllCurve;

public class EllFracSimplifier<K> : TriVarFracSimplifier<K> where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public TriVarPoly<K> EqEll { get; }
    public TriVarPoly<K> DivPol { get; }
    public TriVarPoly<K> SD { get; }

    public EllFracSimplifier(TriVarPoly<K> eqEll, TriVarPoly<K> sd, TriVarPoly<K> dvp)
    {
        (EqEll,SD, DivPol) = (eqEll, sd, dvp);
        // Console.WriteLine(this);
    }

    public override bool IsDivZero(TriVarFrac<K> P)
    {
        if (P.IsZero())
            return true;

        if (DivPol.Degree == 0)
            return false;

        var fdiv = DivPol;
        var arr = P.Num.DecomposeX2().Values.Where(e => !e.IsZero());
        if (arr.All(f => f.Degree > 0 && Ring.FastGCD(f, fdiv).Degree > 0))
            return true;

        return false;
    }

    public override TriVarFrac<K> Apply(TriVarPoly<K> num, TriVarPoly<K> denom)
    {
        if (denom.DegreeOfX2 % 2 == 1)
        {
            num *= SD;
            denom *= SD;
        }

        num = num.Div(EqEll).rem;
        denom = denom.Div(EqEll).rem;
        if (DivPol.Degree > 0)
        {
            num = num.Div(DivPol).rem;
            denom = denom.Div(DivPol).rem;
        }
        else
        {
            var decNum = num.DecomposeX2();
            var decDenom = denom.DecomposeX2();
            var arr = decNum.Values.Concat(decDenom.Values).Where(e => e.Degree > 0).Distinct().ToArray();
            if (arr.Length > 1)
            {
                var gcd = Ring.FastGCD(arr).Monic;
                num = decNum.Select(e => e.Key * (e.Value / gcd)).Aggregate(gcd.Zero, (acc, e) => e + acc);
                denom = decDenom.Select(e => e.Key * (e.Value / gcd)).Aggregate(gcd.Zero, (acc, e) => e + acc);
            }
        }

        return new(this, num, denom);
    }
    
    public override string ToString() => $"EllFracSimplifier[{typeof(K)}][{GetHashCode()}]";
}