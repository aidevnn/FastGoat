using System.Numerics;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.UserGroup.Lattice;

public readonly struct BGVCipher : IElt<BGVCipher>, IRingElt<BGVCipher>, IModuleElt<Rq, BGVCipher>
{
    public Rq A { get; }
    public Rq B { get; }
    public Rq PM { get; }
    public Rational Q { get; }
    public (Rq pm, Rational q) PM_Q => (PM, Q);

    public BGVCipher(Rq a, Rq b, Rq pm, Rational q)
    {
        (A, B, PM, Q) = (a, b, pm, q);
        Hash = HashCode.Combine("bgvcipher", A, B);
    }

    public BGVCipher CoefsMod(Rational q) => (A.CoefsMod(q), B.CoefsMod(q), PM, q);

    public void Show(string name = "")
    {
        if (!string.IsNullOrEmpty(name))
            Console.WriteLine(name);
        Console.WriteLine($"Q:{Q}    PM:{PM}");
        Console.WriteLine($"A:{A}");
        Console.WriteLine($"B:{B}");
    }

    public bool Equals(BGVCipher other) =>
        Q.Equals(other.Q) && PM.Equals(other.PM) && A.Equals(other.A) && B.Equals(other.B);

    public int CompareTo(BGVCipher other)
    {
        return (PM, Q, A, B).CompareTo((other.PM, other.Q, other.A, other.B));
    }

    public int Hash { get; }

    public bool IsZero() => A.CoefsMod(Q).IsZero() && B.CoefsMod(Q).IsZero();

    public Rq KZero => PM.Zero;
    public Rq KOne => PM.One;

    public BGVCipher KMul(Rq k) => new((k * A).ResMod(PM, Q), (k * B).ResMod(PM, Q), PM, Q);

    public BGVCipher Zero => new(PM.Zero, PM.Zero, PM, Q);
    public BGVCipher One => new(PM.One, PM.Zero, PM, Q);

    public BGVCipher Add(BGVCipher e) => new((A + e.A).CoefsMod(Q), (B + e.B).CoefsMod(Q), PM, Q);

    public BGVCipher Sub(BGVCipher e) => new((A - e.A).CoefsMod(Q), (B - e.B).CoefsMod(Q), PM, Q);

    public BGVCipher Opp() => new(A.Opp().CoefsMod(Q), B.Opp().CoefsMod(Q), PM, Q);

    public BGVCipher Mul(BGVCipher e) => new((A * e.A).ResMod(PM, Q), (B * e.B).ResMod(PM, Q), PM, Q);

    public (BGVCipher quo, BGVCipher rem) Div(BGVCipher e)
    {
        throw new NotImplementedException();
    }

    public BGVCipher Mul(int k) => KMul(k * KOne);

    public BGVCipher Pow(int k)
    {
        throw new NotImplementedException();
    }

    public override int GetHashCode() => Hash;

    public override string ToString() => $"(A:{A}, B:{B}) mod {Q} mod {PM}";

    public static BGVCipher operator +(BGVCipher a, BGVCipher b) => a.Add(b);

    public static BGVCipher operator +(int a, BGVCipher b) => b.One.Mul(a).Add(b);

    public static BGVCipher operator +(BGVCipher a, int b) => b + a;

    public static BGVCipher operator -(BGVCipher a) => a.Opp();

    public static BGVCipher operator -(BGVCipher a, BGVCipher b) => a.Sub(b);

    public static BGVCipher operator -(int a, BGVCipher b) => a + b.Opp();

    public static BGVCipher operator -(BGVCipher a, int b) => a + (-b);

    public static BGVCipher operator *(BGVCipher a, BGVCipher b) => a.Mul(b);

    public static BGVCipher operator *(int a, BGVCipher b) => b.Mul(a);

    public static BGVCipher operator *(BGVCipher a, int b) => a.Mul(b);

    public static BGVCipher operator /(BGVCipher a, BGVCipher b)
    {
        throw new NotImplementedException();
    }

    public static BGVCipher operator /(BGVCipher a, int b)
    {
        throw new NotImplementedException();
    }

    public static BGVCipher operator +(BGVCipher a, Rq b) => new((a.A + b).ResMod(a.PM, a.Q), a.B, a.PM, a.Q);

    public static BGVCipher operator +(Rq a, BGVCipher b) => b + a;

    public static BGVCipher operator -(BGVCipher a, Rq b) => (-a) + b;

    public static BGVCipher operator -(Rq a, BGVCipher b) => b + (-a);

    public static BGVCipher operator *(BGVCipher a, Rq b) => a.KMul(b);

    public static BGVCipher operator *(Rq a, BGVCipher b) => b.KMul(a);

    public static BGVCipher operator +(BGVCipher a, Rational b) => a + (b * a.KOne);

    public static BGVCipher operator +(Rational a, BGVCipher b) => b + a;

    public static BGVCipher operator -(BGVCipher a, Rational b) => a + (-b);

    public static BGVCipher operator -(Rational a, BGVCipher b) => a + (-b);

    public static BGVCipher operator *(BGVCipher a, Rational b) => a.KMul(b * a.KOne);

    public static BGVCipher operator *(Rational a, BGVCipher b) => b * a;

    public static implicit operator BGVCipher((Rq a, Rq b, Rq pm, Rational q) e) => new(e.a, e.b, e.pm, e.q);
    public static implicit operator (Rq a, Rq b, Rq pm, Rational q)(BGVCipher e) => (e.A, e.B, e.PM, e.Q);
}