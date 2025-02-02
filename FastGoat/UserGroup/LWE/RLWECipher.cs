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
    public Rational Q { get; }
    public (Rq pm, Rational t, Rational q) PM_T_Q => (PM, T, Q);
    public string Params => $"RLWECipher Q:{Q}    T:{T}    PM:{PM}";

    public RLWECipher(Rq a, Rq b, Rq pm, Rational t, Rational q)
    {
        (A, B, PM, T, Q) = (a, b, pm, t, q);
        Hash = HashCode.Combine("rlwe", A, B);
    }

    public RLWECipher Clone() => new(A.Coefs.ToKPoly(), B.Coefs.ToKPoly(), PM, T, Q);
    public RLWECipher Trunc() => new(A.TruncPoly(), B.TruncPoly(), PM, T, Q);
    public RLWECipher Round() => new(A.RoundPoly(), B.RoundPoly(), PM, T, Q);

    public RLWECipher ModSwitch(Rational qf)
    {
        var f = qf / Q;
        var a = (A * f).ClosestModulusTo(A, T).CoefsModSigned(qf);
        var b = (B * f).ClosestModulusTo(B, T).CoefsModSigned(qf);
        return new(a, b, PM, T, qf);
    }
    
    public RLWECipher ClosestModulusTo(RLWECipher dest)
    {
        var a = A.ClosestModulusTo(dest.A, T).CoefsModSigned(Q);
        var b = B.ClosestModulusTo(dest.B, T).CoefsModSigned(Q);
        return new(a, b, PM, T, Q);
    }
    
    public RLWECipher CoefsModSigned(Rational q0)
    {
        return new RLWECipher(A.CoefsModSigned(q0), B.CoefsModSigned(q0), PM, T, q0);
    }

    public bool Equals(RLWECipher other) =>
        Q.Equals(other.Q) && 
        PM.Equals(other.PM) && A.Equals(other.A) && B.Equals(other.B);

    public int CompareTo(RLWECipher other)
    {
        return (PM_T_Q, A, B).CompareTo((other.PM_T_Q, other.A, other.B));
    }

    public void Show(string name = "")
    {
        if (!string.IsNullOrEmpty(name))
            Console.WriteLine(name);
        Console.WriteLine(Params);
        Console.WriteLine($"A:{A}");
        Console.WriteLine($"B:{B}");
    }

    public Rational KZero => Rational.KZero();
    public Rational KOne => Rational.KOne();
    public int Hash { get; }
    public bool IsZero() => A.IsZero() && B.IsZero();
    public RLWECipher Zero => new(PM.Zero, PM.Zero, PM, T, Q);
    public RLWECipher One => new(PM.One, PM.Zero, PM, T, Q);

    public RLWECipher Add(RLWECipher e) =>
        new((A + e.A).CoefsModSigned(Q), (B + e.B).CoefsModSigned(Q), PM, T, Q);
    public RLWECipher Sub(RLWECipher e) => 
        new((A - e.A).CoefsModSigned(Q), (B - e.B).CoefsModSigned(Q), PM, T, Q);
    public RLWECipher Opp() => new(-A, -B, PM, T, Q);

    public RLWECipher Mul(RLWECipher e) =>
        new((A * e.A).ResModSigned(PM, Q), (B * e.B).ResModSigned(PM, Q), PM, T, Q);

    public (RLWECipher quo, RLWECipher rem) Div(RLWECipher e)
    {
        throw new NotImplementedException();
    }

    public RLWECipher Mul(int k) => new((A * k).CoefsMod(Q), (B * k).CoefsMod(Q), PM, T, Q);

    public RLWECipher Pow(int k)
    {
        throw new NotImplementedException();
    }

    public RLWECipher this[int index] => new(A[index] * PM.One, B[index] * PM.One, PM, T, Q);
    public RLWECipher KMul(Rational k) => new((A * k).CoefsMod(Q), (B * k).CoefsMod(Q), PM, T, Q);
    public override string ToString() => $"(A:{A}, B:{B}) mod {PM} mod {Q} T = {T}";

    public void Deconstruct(out Rq a, out Rq b, out Rq pm, out Rational t, out Rational q)
    {
        (a, b, pm, t, q) = (A, B, PM, T, Q);
    }

    public static RLWECipher operator +(RLWECipher a, RLWECipher b) => a.Add(b);

    public static RLWECipher operator +(int a, RLWECipher b)
    {
        return new RLWECipher((b.A + a).CoefsMod(b.Q), b.B, b.PM, b.T, b.Q);
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
        new((a.A + b).CoefsMod(a.Q), a.B, a.PM, a.T, a.Q);

    public static RLWECipher operator +(Rational a, RLWECipher b) => b + a;

    public static RLWECipher operator -(RLWECipher a, Rational b) => a + (-b);

    public static RLWECipher operator -(Rational a, RLWECipher b) => (-b) + a;

    public static RLWECipher operator *(RLWECipher a, Rational b) => a.KMul(b);

    public static RLWECipher operator *(Rational a, RLWECipher b) => b * a;

    public static RLWECipher operator +(RLWECipher a, Rq b) =>
        new((a.A + b).CoefsMod(a.Q), a.B, a.PM, a.T, a.Q);

    public static RLWECipher operator +(Rq a, RLWECipher b) => b + a;

    public static RLWECipher operator -(RLWECipher a, Rq b) => a + (-b);

    public static RLWECipher operator -(Rq a, RLWECipher b) => (-b) + a;

    public static RLWECipher operator *(RLWECipher a, Rq b)
    {
        var ca = (a.A * b).ResModSigned(a.PM, a.Q);
        var cb = (a.B * b).ResModSigned(a.PM, a.Q);
        return new RLWECipher(ca, cb, a.PM, a.T, a.Q);
    }

    public static RLWECipher operator *(Rq a, RLWECipher b) => b * a;

    public static implicit operator RLWECipher((Rq a, Rq b, Rq pm, Rational t, Rational q) e)
        => new(e.a, e.b, e.pm, e.t, e.q);
}