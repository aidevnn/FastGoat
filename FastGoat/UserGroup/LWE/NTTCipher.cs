using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Padic;

namespace FastGoat.UserGroup.LWE;

public readonly struct NTTCipher : IModuleElt<KMatrix<ZnBigInt>, NTTCipher>, IElt<NTTCipher>, IRingElt<NTTCipher>
{
    public KMatrix<ZnBigInt> A { get; }

    public KMatrix<ZnBigInt> B { get; }
    public NTTInfos NttInfos { get; }

    public NTTCipher(KMatrix<ZnBigInt> a, KMatrix<ZnBigInt> b, NTTInfos nttInfos)
    {
        (A, B, NttInfos) = (a, b, nttInfos);
        Hash = (A.Hash, B.Hash).GetHashCode();
    }

    public bool Equals(NTTCipher other) => (A, B).Equals((other.A, other.B));

    public int CompareTo(NTTCipher other) => (A, B).CompareTo((other.A, other.B));

    public KMatrix<ZnBigInt> KZero => A.Zero;
    public KMatrix<ZnBigInt> KOne => NttInfos.ntt.GetCol(0);

    public NTTCipher KMul(KMatrix<ZnBigInt> k)
    {
        return new(RLWE.MulNTT(A, k), RLWE.MulNTT(B, k), NttInfos);
    }

    public static NTTCipher operator +(NTTCipher a, KMatrix<ZnBigInt> b)
    {
        return new(a.A + b, a.B, a.NttInfos);
    }

    public static NTTCipher operator +(KMatrix<ZnBigInt> a, NTTCipher b) => b + a;

    public static NTTCipher operator -(NTTCipher a, KMatrix<ZnBigInt> b) => a + (-b);

    public static NTTCipher operator -(KMatrix<ZnBigInt> a, NTTCipher b) => a + (-b);

    public static NTTCipher operator *(NTTCipher a, KMatrix<ZnBigInt> b) => a.KMul(b);

    public static NTTCipher operator *(KMatrix<ZnBigInt> a, NTTCipher b) => b * a;

    public int Hash { get; }
    public bool IsZero() => A.IsZero() && B.IsZero();

    public NTTCipher Clone()
    {
        var nttInfos = new NTTInfos(NttInfos.n, NttInfos.w, NttInfos.ntt.Clone, NttInfos.wPows.ToArray(), NttInfos.t,
            NttInfos.intt.Clone, NttInfos.iwPows.ToArray());
        var a = A.Clone;
        var b = B.Clone;
        return new(a, b, nttInfos);
    }

    public NTTCipher Zero => new(A.Zero, B.Zero, NttInfos);
    public NTTCipher One => new(KOne, B.Zero, NttInfos);
    public NTTCipher ESK => new(A.Zero, -KOne, NttInfos);
    public NTTCipher Add(NTTCipher e) => new(A + e.A, B + e.B, NttInfos);

    public NTTCipher Sub(NTTCipher e) => new(A - e.A, B - e.B, NttInfos);

    public NTTCipher Opp() => new(-A, -B, NttInfos);

    public NTTCipher Mul(NTTCipher e)
    {
        throw new NotImplementedException();
    }

    public (NTTCipher quo, NTTCipher rem) Div(NTTCipher e)
    {
        throw new NotImplementedException();
    }

    public NTTCipher Mul(int k) => KMul(k * KOne);
    public NTTCipher Mul(ZnBigInt k) => KMul(k * KOne);

    public NTTCipher Pow(int k)
    {
        throw new NotImplementedException();
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        return $"A:{A.T} B:{B.T} mod {KOne.KOne.Mod}";
    }

    public static NTTCipher operator +(NTTCipher a, NTTCipher b) => a.Add(b);

    public static NTTCipher operator +(int a, NTTCipher b) => a * b.KOne + b;

    public static NTTCipher operator +(NTTCipher a, int b) => b + a;

    public static NTTCipher operator -(NTTCipher a) => a.Opp();

    public static NTTCipher operator -(NTTCipher a, NTTCipher b) => a.Sub(b);

    public static NTTCipher operator -(int a, NTTCipher b) => a + (-b);

    public static NTTCipher operator -(NTTCipher a, int b) => a + (-b);

    public static NTTCipher operator *(NTTCipher a, NTTCipher b)
    {
        throw new NotImplementedException();
    }

    public static NTTCipher operator *(int a, NTTCipher b) => b.Mul(a);

    public static NTTCipher operator *(NTTCipher a, int b) => a.Mul(b);

    public static NTTCipher operator *(ZnBigInt a, NTTCipher b) => b.Mul(a);

    public static NTTCipher operator *(NTTCipher a, ZnBigInt b) => a.Mul(b);

    public static NTTCipher operator /(NTTCipher a, NTTCipher b)
    {
        throw new NotImplementedException();
    }

    public static NTTCipher operator /(NTTCipher a, int b) => a * (a.KOne * (a.A.KOne * b).Inv());

    public static implicit operator NTTCipher((KMatrix<ZnBigInt> a, KMatrix<ZnBigInt> b, NTTInfos tables) e) =>
        new(e.a, e.b, e.tables);
}