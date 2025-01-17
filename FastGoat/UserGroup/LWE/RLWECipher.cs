using FastGoat.Structures;
using FastGoat.UserGroup.Integers;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.UserGroup.LWE;

public struct RLWECipher : IModuleElt<Rational, RLWECipher>, IElt<RLWECipher>, IRingElt<RLWECipher>
{
    public Rq A { get; }
    public Rq B { get; }
    public Rq PM { get; }
    public Rational Q { get; }
    public (Rq pm, Rational q) PM_Q => (PM, Q);

    public RLWECipher(Rq a, Rq b, Rq pm, Rational q)
    {
        (A, B, PM, Q) = (a, b, pm, q);
        Hash = HashCode.Combine("rlwe", A, B);
    }

    public RLWECipher CoefsMod(Rational p) => new(A.ResMod(PM, p), B.ResMod(PM, p), PM, p);
    public RLWECipher CoefsModSigned(Rational t) => new(A.CoefsModSigned(t), B.CoefsModSigned(t), PM, t);
    public RLWECipher CoefsModSigned() => CoefsModSigned(Q);
    public RLWECipher Trunc() => new(A.TruncPoly(), B.TruncPoly(), PM, Q);
    public RLWECipher Round() => new(A.RoundPoly(), B.RoundPoly(), PM, Q);

    public bool Equals(RLWECipher other) =>
        Q.Equals(other.Q) && PM.Equals(other.PM) && A.Equals(other.A) && B.Equals(other.B);

    public int CompareTo(RLWECipher other)
    {
        return (PM, Q, A, B).CompareTo((other.PM, other.Q, other.A, other.B));
    }

    public void Show(string name = "")
    {
        if (!string.IsNullOrEmpty(name))
            Console.WriteLine(name);
        Console.WriteLine($"Q:{Q}    PM:{PM}");
        Console.WriteLine($"A:{A}");
        Console.WriteLine($"B:{B}");
    }

    public Rational KZero => Rational.KZero();
    public Rational KOne => Rational.KOne();
    public int Hash { get; }
    public bool IsZero() => A.IsZero() && B.IsZero();

    public RLWECipher Zero => new(PM.Zero, PM.Zero, PM, Q);
    public RLWECipher One => new(PM.One, PM.Zero, PM, Q);
    public RLWECipher Add(RLWECipher e) => new((A + e.A).CoefsMod(Q), (B + e.B).CoefsMod(Q), PM, Q);

    public RLWECipher Sub(RLWECipher e) => new((A - e.A).CoefsMod(Q), (B - e.B).CoefsMod(Q), PM, Q);

    public RLWECipher Opp() => new(-A, -B, PM, Q);

    public RLWECipher Mul(RLWECipher e) => new((A * e.A).CoefsMod(Q), (B * e.B).CoefsMod(Q), PM, Q);

    public (RLWECipher quo, RLWECipher rem) Div(RLWECipher e)
    {
        throw new NotImplementedException();
    }

    public RLWECipher Mul(int k) => new((A * k).CoefsMod(Q), (B * k).CoefsMod(Q), PM, Q);

    public RLWECipher Pow(int k)
    {
        throw new NotImplementedException();
    }

    public RLWECipher this[int index] => new(A[index] * PM.One, B[index] * PM.One, PM, Q);

    public RLWECipher KMul(Rational k) => new((A * k).CoefsMod(Q), (B * k).CoefsMod(Q), PM, Q);

    public override string ToString() => $"(A:{A}, B:{B}) mod {Q} mod {PM}";

    public void Deconstruct(out Rq a, out Rq b, out Rq pm, out Rational q)
    {
        (a, b, pm, q) = (A, B, PM, Q);
    }

    public static RLWECipher operator +(RLWECipher a, RLWECipher b) => a.Add(b);

    public static RLWECipher operator +(int a, RLWECipher b) => new((b.A + a).CoefsMod(b.Q), b.B, b.PM, b.Q);

    public static RLWECipher operator +(RLWECipher a, int b) => b + a;

    public static RLWECipher operator -(RLWECipher a) => a.Opp();

    public static RLWECipher operator -(RLWECipher a, RLWECipher b) => a.Sub(b);

    public static RLWECipher operator -(int a, RLWECipher b) => a + (-b);

    public static RLWECipher operator -(RLWECipher a, int b) => a + (-b);

    public static RLWECipher operator *(RLWECipher a, RLWECipher b) => a.Mul(b);

    public static RLWECipher operator *(int a, RLWECipher b) => b.Mul(a);

    public static RLWECipher operator *(RLWECipher a, int b) => a.Mul(b);

    public static RLWECipher operator /(RLWECipher a, RLWECipher b)
    {
        throw new NotImplementedException();
    }

    public static RLWECipher operator /(RLWECipher a, int b)
    {
        throw new NotImplementedException();
    }

    public static RLWECipher operator +(RLWECipher a, Rational b) =>
        new((a.A + b).CoefsMod(a.Q), a.B.CoefsMod(a.Q), a.PM, a.Q);

    public static RLWECipher operator +(Rational a, RLWECipher b) => b + a;

    public static RLWECipher operator -(RLWECipher a, Rational b) => a + (-b);

    public static RLWECipher operator -(Rational a, RLWECipher b) => (-b) + a;

    public static RLWECipher operator *(RLWECipher a, Rational b) => a.KMul(b);

    public static RLWECipher operator *(Rational a, RLWECipher b) => b * a;

    public static RLWECipher operator +(RLWECipher a, Rq b) =>
        new((a.A + b).CoefsMod(a.Q), a.B.CoefsMod(a.Q), a.PM, a.Q);

    public static RLWECipher operator +(Rq a, RLWECipher b) => b + a;

    public static RLWECipher operator -(RLWECipher a, Rq b) => a + (-b);

    public static RLWECipher operator -(Rq a, RLWECipher b) => (-b) + a;

    public static RLWECipher operator *(RLWECipher a, Rq b) =>
        new((a.A * b).CoefsMod(a.Q), (a.B * b).CoefsMod(a.Q), a.PM, a.Q);

    public static RLWECipher operator *(Rq a, RLWECipher b) => b * a;

    public static implicit operator RLWECipher((Rq a, Rq b, Rq pm, Rational q) e) => new(e.a, e.b, e.pm, e.q);
}