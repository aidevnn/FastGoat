using FastGoat.Structures;
using FastGoat.UserGroup.Integers;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.UserGroup.LWE;

public struct RLWECipher : IModuleElt<Rational, RLWECipher>, IElt<RLWECipher>, IRingElt<RLWECipher>
{
    public Rq A { get; }
    public Rq B { get; }
    public Rq PM { get; }
    public Rational T { get; }
    public Rational Q0 { get; }
    public Rational Q1 { get; }
    public Rational Q2 { get; }
    public (Rq pm, Rational t, Rational q0, Rational q1, Rational q2) PM_T_Q => (PM, T, Q0, Q1, Q2);

    public RLWECipher(Rq a, Rq b, Rq pm, Rational t, Rational q0, Rational q1, Rational q2)
    {
        (A, B, PM, T, Q0, Q1, Q2) = (a, b, pm, t, q0, q1, q2);
        Hash = HashCode.Combine("rlwe", A, B);
    }

    public RLWECipher Clone() => new(A.Coefs.ToKPoly(), B.Coefs.ToKPoly(), PM, T, Q0, Q1, Q2);
    public RLWECipher Trunc() => new(A.TruncPoly(), B.TruncPoly(), PM, T, Q0, Q1, Q2);
    public RLWECipher Round() => new(A.RoundPoly(), B.RoundPoly(), PM, T, Q0, Q1, Q2);

    public RLWECipher ModSwitch(Rational qi, Rational qf)
    {
        var f = qf / qi;
        var a = (A * f).ClosestModulusTo(A, T).CoefsModSigned(qf);
        var b = (B * f).ClosestModulusTo(B, T).CoefsModSigned(qf);
        return new(a, b, PM, T, Q0, Q1, Q2);
    }

    public bool Equals(RLWECipher other) =>
        Q0.Equals(other.Q0) && Q1.Equals(other.Q1) && Q2.Equals(other.Q2) &&
        PM.Equals(other.PM) && A.Equals(other.A) && B.Equals(other.B);

    public int CompareTo(RLWECipher other)
    {
        return (PM_T_Q, A, B).CompareTo((other.PM_T_Q, other.A, other.B));
    }

    public void Show(string name = "")
    {
        if (!string.IsNullOrEmpty(name))
            Console.WriteLine(name);
        Console.WriteLine($"Q2:{Q2}    Q1:{Q1}    Q0:{Q0}    T:{T}    PM:{PM}");
        Console.WriteLine($"A:{A}");
        Console.WriteLine($"B:{B}");
    }

    public Rational KZero => Rational.KZero();
    public Rational KOne => Rational.KOne();
    public int Hash { get; }
    public bool IsZero() => A.IsZero() && B.IsZero();
    public RLWECipher Zero => new(PM.Zero, PM.Zero, PM, T, Q0, Q1, Q2);
    public RLWECipher One => new(PM.One, PM.Zero, PM, T, Q0, Q1, Q2);

    public RLWECipher Add(RLWECipher e) =>
        new((A + e.A).CoefsModSigned(Q1), (B + e.B).CoefsModSigned(Q1), PM, T, Q0, Q1, Q2);
    public RLWECipher Sub(RLWECipher e) => 
        new((A - e.A).CoefsModSigned(Q1), (B - e.B).CoefsModSigned(Q1), PM, T, Q0, Q1, Q2);
    public RLWECipher Opp() => new(-A, -B, PM, T, Q0, Q1, Q2);

    public RLWECipher Mul(RLWECipher e) =>
        new((A * e.A).ResModSigned(PM, Q1), (B * e.B).ResModSigned(PM, Q1), PM, T, Q0, Q1, Q2);

    public (RLWECipher quo, RLWECipher rem) Div(RLWECipher e)
    {
        throw new NotImplementedException();
    }

    public RLWECipher Mul(int k) => new((A * k).CoefsMod(Q1), (B * k).CoefsMod(Q1), PM, T, Q0, Q1, Q2);

    public RLWECipher Pow(int k)
    {
        throw new NotImplementedException();
    }

    public RLWECipher this[int index] => new(A[index] * PM.One, B[index] * PM.One, PM, T, Q0, Q1, Q2);
    public RLWECipher KMul(Rational k) => new((A * k).CoefsMod(Q1), (B * k).CoefsMod(Q1), PM, T, Q0, Q1, Q2);
    public override string ToString() => $"(A:{A}, B:{B}) mod {PM} mod {Q1} T = {T}";

    public void Deconstruct(out Rq a, out Rq b, out Rq pm, out Rational t, out Rational q0, out Rational q1,
        out Rational q2)
    {
        (a, b, pm, t, q0, q1, q2) = (A, B, PM, T, Q0, Q1, Q2);
    }

    public static RLWECipher operator +(RLWECipher a, RLWECipher b) => a.Add(b);

    public static RLWECipher operator +(int a, RLWECipher b)
    {
        return new RLWECipher((b.A + a).CoefsMod(b.Q0), b.B, b.PM, b.T, b.Q0, b.Q1, b.Q2);
    }

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
        new((a.A + b).CoefsMod(a.Q1), a.B, a.PM, a.T, a.Q0, a.Q1, a.Q2);

    public static RLWECipher operator +(Rational a, RLWECipher b) => b + a;

    public static RLWECipher operator -(RLWECipher a, Rational b) => a + (-b);

    public static RLWECipher operator -(Rational a, RLWECipher b) => (-b) + a;

    public static RLWECipher operator *(RLWECipher a, Rational b) => a.KMul(b);

    public static RLWECipher operator *(Rational a, RLWECipher b) => b * a;

    public static RLWECipher operator +(RLWECipher a, Rq b) =>
        new((a.A + b).CoefsMod(a.Q1), a.B, a.PM, a.T, a.Q0, a.Q1, a.Q2);

    public static RLWECipher operator +(Rq a, RLWECipher b) => b + a;

    public static RLWECipher operator -(RLWECipher a, Rq b) => a + (-b);

    public static RLWECipher operator -(Rq a, RLWECipher b) => (-b) + a;

    public static RLWECipher operator *(RLWECipher a, Rq b)
    {
        var ca = (a.A * b).ResModSigned(a.PM, a.Q1);
        var cb = (a.B * b).ResModSigned(a.PM, a.Q1);
        return new RLWECipher(ca, cb, a.PM, a.T, a.Q0, a.Q1, a.Q2);
    }

    public static RLWECipher operator *(Rq a, RLWECipher b) => b * a;

    public static implicit operator RLWECipher((Rq a, Rq b, Rq pm, Rational t, Rational q0, Rational q1, Rational q2) e)
        => new(e.a, e.b, e.pm, e.t, e.q0, e.q1, e.q2);
}