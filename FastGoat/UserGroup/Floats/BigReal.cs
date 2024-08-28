using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Floats;

public enum Rounding
{
    SciForm,
    FixForm
}

public readonly struct BigReal : IElt<BigReal>, IRingElt<BigReal>, IFieldElt<BigReal>, IVsElt<Rational, BigReal>, IFloatElt<BigReal>
{
    public static DigitsForm Display = DigitsForm.Default;

    public enum DigitsForm
    {
        Default,
        SciForm
    }

    public static BigReal BrZero(int o) => o > 0 ? new(0, 0, o) : throw new();
    public static BigReal BrOne(int o) => o > 0 ? new(1, 0, o) : throw new();
    public static BigReal BrPow10n(int n, int o) => o > 0 ? new(1, n, o) : throw new();
    public static BigReal Pi(int O) => FromBigInteger(BigInteger.Parse(PiStr.Take(O).Glue()), O).Mul10PowN(-O + 1);
    public static BigReal E(int O) => FromBigInteger(BigInteger.Parse(EStr.Take(O).Glue()), O).Mul10PowN(-O + 1);
    public BigInteger K { get; }
    public int NbDigits { get; }
    public int V { get; }
    public int O { get; }

    public (BigInteger K, int V, int NbDigits, int O) Details => (K, V, NbDigits, O);

    private BigReal(BigInteger k, int v, int o)
    {
        K = k;
        V = v;
        if (k.IsZero)
            V = 0;

        NbDigits = Length(K);
        O = o;
        Hash = (K, V, O).GetHashCode();
    }

    public bool Equals(BigReal other) => Sub(other).IsZero();

    public int CompareTo(BigReal other)
    {
        var sub = Sub(other);
        return sub.IsZero() ? 0 : sub.K.Sign;
    }

    public int Hash { get; }
    public bool IsZero() => K == 0 || V < -O + 2;
    public bool IsZero3d() => K == 0 || V < -O + 3;
    public bool IsZero4d() => K == 0 || V < -O + 4;

    public BigReal Zero => new(0, 0, O);
    public BigReal One => new(BigInteger.Pow(10, O - 1), 0, O);
    public int P => 0;
    public bool Invertible() => !IsZero();

    public BigReal Pow10(int n) => new(BigInteger.Pow(10, O - 1), n, O);
    public BigReal Mul10PowN(int N) => new(K, V + N, O);

    public BigReal Add(BigReal e)
    {
        if (O < e.O)
            return ToBigReal(e.O).Add(e);

        var d = V - e.V;
        var ad = Int32.Abs(d);
        if (ad > O)
        {
            if (d < 0)
                return e;
            else
                return this;
        }
        else
        {
            if (d >= 0)
            {
                if (NbDigits >= e.NbDigits)
                {
                    var shift0 = BigInteger.Pow(10, d);
                    var shift1 = BigInteger.Pow(10, NbDigits - e.NbDigits);
                    var k0 = e.K * shift1 + shift0 * K;
                    var k1 = Clamp(k0, O);
                    var a = d + NbDigits - Length(k0);
                    var r = new BigReal(k1, V - a, O);
                    return r;
                }
                else
                {
                    var shift0 = BigInteger.Pow(10, d);
                    var shift1 = BigInteger.Pow(10, e.NbDigits - NbDigits);
                    var k0 = e.K + shift1 * shift0 * K;
                    var k1 = Clamp(k0, O);
                    var a = d + e.NbDigits - Length(k0);
                    var r = new BigReal(k1, V - a, O);
                    return r;
                }
            }
            else
                return e.Add(this);
        }
    }

    public BigReal Opp() => new(-K, V, O);

    public BigReal Sub(BigReal e) => Add(e.Opp());

    public BigReal Mul(BigReal e)
    {
        if (IsZero() || e.IsZero())
            return Zero;

        if (O < e.O)
            return ToBigReal(e.O).Mul(e);

        var k0 = K * e.K;
        var n0 = V + e.V;
        var a = Length(k0) - NbDigits - e.NbDigits + 1;
        var k1 = Clamp(k0, O);
        return new(k1, n0 + a, O);
    }

    public (BigReal quo, BigReal rem) Div(BigReal e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();

        if (O < e.O)
            return ToBigReal(e.O).Div(e);

        var k0 = (K * BigInteger.Pow(10, O + 4)) / e.K;
        var d = Length(k0) - (O + NbDigits - e.NbDigits + 5);
        var k1 = Clamp(k0, O);
        var n1 = V - e.V;
        return (new(k1, n1 + d, O), Zero);
    }

    public BigReal Mul(int k)
    {
        var e = FromBigInteger(k, O);
        return Mul(e);
    }

    public BigReal Inv() => One.Div(this).quo;

    public BigReal Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        var br = One;
        for (int i = 0; i < k; i++)
            br *= this;

        return br;
    }

    public override int GetHashCode() => Hash;

    public BigReal ToBigReal(int o)
    {
        if (V < 0 && -V > o)
            return new(0, 0, o);

        var k0 = Clamp(K, o);
        return new(k0, V, o);
    }

    public BigReal Round0 => Round(this);

    public double ToDouble => double.Parse(ToSciForm());
    
    public decimal ToDecimal => decimal.Parse(ToFixForm());

    public BigReal Absolute => new(BigInteger.Abs(K), V, O);
    public BigReal Sqrt() => Sqrt(this);

    public Rational ToRational
    {
        get
        {
            if (V - NbDigits + 1 >= 0)
            {
                var num = K * BigInteger.Pow(10, V - NbDigits + 1);
                return new(num);
            }
            else
            {
                var num = K;
                var denom = BigInteger.Pow(10, NbDigits - V - 1);
                return new(num, denom);
            }
        }
    }

    public int Sign => K.Sign;
    public Dcml ToDcml => new(ToDecimal);
    public Dble ToDble => new(ToDouble);
    
    public static BigReal FromFixedPrecision<T>(T e, int O) 
        where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>, IFloatElt<T>, IFixedPrecisionElt<T>
    {
        if (e is Dble e0)
            return e0.ToBigReal(O);
        if (e is Dcml e1)
            return e1.ToBigReal(O);

        throw new ArgumentException();
    }

    public string ToSciForm()
    {
        var v = V < 0 ? $"-{-V:000}" : $"+{V:000}";
        var s0 = K < 0 ? $"{-K}" : $"{K}";
        var s1 = $"{s0[0]}.{s0.Skip(1).Concat(Enumerable.Repeat('0', O)).Take(O - 1).Glue()}";
        return K < 0 ? $"-{s1}E{v}" : $"{s1}E{v}";
    }

    public string ToFixForm()
    {
        if (-V > O)
            return ToSciForm();

        var k0 = K.Sign == 1 ? this : Opp();
        if (V == 0)
            return ToSciForm().Split("E")[0];
        if (V > 0)
        {
            var s0 = k0.ToSciForm().Split("E")[0].Replace(".", "");
            s0 += Enumerable.Repeat('0', V).Glue();
            s0 = s0.InsertAt(V + 1, '.').Glue();
            s0 = s0.TrimEnd('0').TrimEnd('.');
            return K.Sign == 1 ? s0 : $"-{s0}";
        }
        else
        {
            var s0 = k0.ToSciForm().Split("E")[0].Replace(".", "");
            s0 = Enumerable.Repeat('0', -V).Glue() + s0;
            s0 = s0.InsertAt(1, '.').Glue();
            s0 = s0.TrimEnd('0').TrimEnd('.');
            return K.Sign == 1 ? s0 : $"-{s0}";
        }
    }

    public override string ToString()
    {
        if (Display == DigitsForm.Default)
            return V < 200 ? $"{ToDouble}" : $"{ToBigReal(16).ToSciForm()}";

        if (IsZero())
            return Enumerable.Repeat("0", O - 1).Prepend("0.").Append("E+000").Glue();

        return ToSciForm();
    }

    public static BigReal Min(BigReal a, BigReal b) => a.CompareTo(b) <= 0 ? a : b;
    public static BigReal Max(BigReal a, BigReal b) => a.CompareTo(b) >= 0 ? a : b;

    public static BigReal Sqrt(BigReal r) => NthRoot(r, 2);
    public static BigReal NthRoot(BigInteger r, int n, int o) => NthRoot(new BigReal(r, 0, o), n);

    public static BigReal NthRoot(BigReal r, int n)
    {
        if (n == 0)
            return r.One;

        if (n < 0)
            return NthRoot(r, -n).Inv();

        if (n % 2 == 0 && r.Sign == -1)
            throw new("Even NthRoot must has positive argument");

        if (n == 1)
            return r;

        BigReal ai;
        var aj = 1 + r / n;
        do
        {
            ai = aj;
            var aiPow = ai.Pow(n - 1);
            var num = aiPow * ai - r;
            var denom = n * aiPow;
            aj = ai - num / denom; // Newton iteration
        } while (!(ai - aj).IsZero());

        return aj;
    }
    public BigReal RoundEven => Round0;

    // TODO Fix
    public static BigReal Round(BigReal br, int d = 0, MidpointRounding mid = MidpointRounding.ToEven,
        Rounding form = Rounding.FixForm)
    {
        if (d < 0)
            throw new("Nb digits to round must be positive");

        if (form == Rounding.SciForm && d > br.NbDigits)
            return br;

        if (form == Rounding.FixForm && br.V - 1 >= br.NbDigits)
            return br;

        if (form == Rounding.FixForm && -d > br.V + 1)
            return br.Zero;

        var d0 = form == Rounding.SciForm ? d + 1 : d + 1 + br.V;
        var d1 = form == Rounding.SciForm ? br.V - d : -d;
        var k0 = d0 == 0 ? BigInteger.Zero : Clamp(br.K, d0);
        var br0 = new BigReal(k0, br.V, br.O);
        var den = new BigReal(1, d1, br.O);
        var r = br - br0;
        var rs = r.K.Sign;
        var r0 = r * rs * 2;
        var comp = r0.CompareTo(den);

        if (comp != 0)
        {
            if (comp == -1)
                return br0;
            else
                return br0 + new BigReal(rs, d1, br.O);
        }
        else
        {
            if (mid == MidpointRounding.ToEven)
            {
                if (BigInteger.IsEvenInteger(br0.K))
                    return br0;
                else
                    return br0 + new BigReal(rs, d1, br.O);
            }
            else if (mid == MidpointRounding.ToNegativeInfinity)
            {
                if (rs == 1)
                    return br0;
                else
                    return br0 + new BigReal(rs, d1, br.O);
            }
            else if (mid == MidpointRounding.ToPositiveInfinity)
            {
                if (rs == 1)
                    return br0 + new BigReal(rs, d1, br.O);
                else
                    return br0;
            }
            else if (mid == MidpointRounding.ToZero)
                return br0;
            else
                return br0 + new BigReal(rs, d1, br.O);
        }
    }

    public static BigReal NormN(int n, BigReal[] br)
    {
        if (n < 1)
            throw new($"n={n} must be >= 1");

        var sum = br.Aggregate(br[0].Zero, (acc, b0) => acc + b0.Absolute.Pow(n));
        return NthRoot(sum, n);
    }

    public static int Length(BigInteger k)
    {
        var k0 = BigInteger.Abs(k);

        if (k0 < 10) return 1;
        if (k0 < 100) return 2;
        if (k0 < 1000) return 3;
        if (k0 < 10000) return 4;
        if (k0 < 100000) return 5;
        if (k0 < 1000000) return 6;
        if (k0 < 10000000) return 7;
        if (k0 < 100000000) return 8;
        if (k0 < 1000000000) return 9;
        if (k0 < 10000000000) return 10;
        if (k0 < 100000000000) return 11;
        if (k0 < 1000000000000) return 12;
        if (k0 < 10000000000000) return 13;
        if (k0 < 100000000000000) return 14;
        if (k0 < 1000000000000000) return 15;
        if (k0 < 10000000000000000) return 16;
        if (k0 < 100000000000000000) return 17;
        if (k0 < 1000000000000000000) return 18;
        if (k0 < 10000000000000000000) return 19;

        return k0.ToString().Length;
    }

    public static BigInteger Clamp(BigInteger k, int O)
    {
        if (O < 1)
            throw new($"Nb digits >= 1 but is {O}");

        if (k > 0)
        {
            var length = Length(k);
            if (length <= O)
                return k;

            var n1 = length - O;
            var k1 = k / BigInteger.Pow(10, n1);
            return k1;
        }

        if (k < 0)
        {
            var length = Length(k);
            if (length <= O)
                return k;

            var n1 = length - O;
            var k1 = k / BigInteger.Pow(10, n1);
            return k1;
        }

        return 0;
    }

    public static BigReal Pow10(int n, int o) => o > 0 ? new(1, n, o) : throw new();

    public static BigReal FromDouble(double d, int o)
    {
        var s = $"{d:E17}".Split('E');
        var s10 = s[0].Replace(".", "");
        var k0 = long.Parse(s10);
        var sgn = s[1][0] == '+' ? 1 : -1;
        var exp = int.Parse(s[1].Skip(1).Glue()) * sgn;
        var k1 = Clamp(k0, o);
        return new(k1, exp, o);
    }

    public static BigReal FromDecimal(decimal d, int o)
    {
        var s = $"{d:E28}".Split('E');
        var s10 = s[0].Replace(".", "");
        var k0 = BigInteger.Parse(s10);
        var sgn = s[1][0] == '+' ? 1 : -1;
        var exp = int.Parse(s[1].Skip(1).Glue()) * sgn;
        var k1 = Clamp(k0, o);
        return new(k1, exp, o);
    }

    public static BigReal FromString(string s, int o)
    {
        s = s.Replace('E', 'e').Trim();
        var s1 = s.Split('e');
        var lt = s1[0].Split('.');
        if (lt.Length > 2)
            throw new($"Invalid format1 for {s}, scientific format ex:1.2345E+004");
        else if (lt.Length == 1)
        {
            var k = BigInteger.Parse(lt[0]);
            var k1 = Clamp(k, o);
            if (k1.IsZero)
                return BrZero(o);

            var length = Length(k1);
            var v = s1.Length == 1 ? 0 : int.Parse(s1[1]);
            return new(k1, v + length - 1, o);
        }
        else
        {
            var (lt0, lt1) = (lt[0].TrimStart(), lt[1].TrimEnd());
            var k0 = BigInteger.Parse(lt0);
            var k1 = BigInteger.Parse(lt1);
            var l0 = s[0] == '-' ? lt0.Length - 1 : lt0.Length;
            var l1 = lt1.Length;
            var sgn = k0.IsZero ? (s[0] == '-' ? -1 : 1) : k0.Sign;
            var k = k1 * sgn + k0 * BigInteger.Pow(10, l1);
            var k2 = Clamp(k, o);
            if (k2.IsZero)
                return BrZero(o);

            var v = s1.Length == 1 ? 0 : int.Parse(s1[1]);
            if (k0.IsZero)
            {
                var t0 = lt1.Length - lt1.TrimStart('0').Length;
                v -= t0 + 1;
            }

            return new(k2, v + l0 - 1, o);
        }
    }

    public static BigReal FromBigInteger(BigInteger r, int o)
    {
        var d = Length(r) - 1;
        var k = Clamp(r, o);
        return new(k, d, o);
    }

    public static BigReal FromBigIntegerAndExponent(BigInteger r, int v, int o)
    {
        var k = Clamp(r, o);
        return new(k, v, o);
    }

    public static BigReal FromRational(Rational r, int o)
    {
        var num = FromBigInteger(r.Num, o);
        var denom = FromBigInteger(r.Denom, o);
        return num / denom;
    }

    public static BigReal operator +(BigReal a, BigReal b) => a.Add(b);

    public static BigReal operator +(int a, BigReal b) => b.One.Mul(a).Add(b);

    public static BigReal operator +(BigReal a, int b) => a.Add(a.One.Mul(b));

    public static BigReal operator -(BigReal a) => a.Opp();

    public static BigReal operator -(BigReal a, BigReal b) => a.Sub(b);

    public static BigReal operator -(int a, BigReal b) => b.One.Mul(a).Sub(b);

    public static BigReal operator -(BigReal a, int b) => a.Sub(a.One.Mul(b));

    public static BigReal operator *(BigReal a, BigReal b) => a.Mul(b);

    public static BigReal operator *(int a, BigReal b) => b.Mul(a);

    public static BigReal operator *(BigReal a, int b) => a.Mul(b);

    public static BigReal operator /(BigReal a, BigReal b) => a.Div(b).quo;

    public static BigReal operator /(BigReal a, int b) => a.Div(a.One.Mul(b)).quo;

    public static BigReal operator /(int a, BigReal b) => b.One.Mul(a).Div(b).quo;

    private const int DigitsDouble = 17;

    public static double Abs(BigReal t) => Double.Abs(t.ToDouble);

    public static bool IsValuedField => true;
    public Rational KZero => Rational.KZero();
    public Rational KOne => Rational.KOne();
    public BigReal KMul(Rational k) => Mul(FromRational(k, O));

    public static BigReal operator +(BigReal a, Rational b) => a.Add(FromRational(b, a.O));

    public static BigReal operator +(Rational a, BigReal b) => FromRational(a, b.O).Add(b);

    public static BigReal operator -(BigReal a, Rational b) => a.Sub(FromRational(b, a.O));

    public static BigReal operator -(Rational a, BigReal b) => FromRational(a, b.O).Sub(b);

    public static BigReal operator *(BigReal a, Rational b) => a.KMul(b);

    public static BigReal operator *(Rational a, BigReal b) => b.KMul(a);

    public static BigReal operator /(BigReal a, Rational b) => a.KMul(b.Inv());

    public static bool operator <(BigReal a, BigReal b) => a.CompareTo(b) == -1;
    public static bool operator >(BigReal a, BigReal b) => a.CompareTo(b) == 1;

    // mathematica cloud
    private const string PiStr =
        "3141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521384146951941511609433057270365759591953092186117381932611793105118548074462379962749567351885752724891227938183011949129833673362440656643086021394946395224737190702179860943702770539217176293176752384674818467669405132000568127145263560827785771342757789609173637178721468440901224953430146549585371050792279689258923542019956112129021960864034418159813629774771309960518707211349999998372978049951059731732816096318595024459455346908302642522308253344685035261931188171010003137838752886587533208381420617177669147303598253490428755468731159562863882353787593751957781857780532171226806613001927876611195909216420199";

    private const string EStr =
        "2718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260595630738132328627943490763233829880753195251019011573834187930702154089149934884167509244761460668082264800168477411853742345442437107539077744992069551702761838606261331384583000752044933826560297606737113200709328709127443747047230696977209310141692836819025515108657463772111252389784425056953696770785449969967946864454905987931636889230098793127736178215424999229576351482208269895193668033182528869398496465105820939239829488793320362509443117301238197068416140397019837679320683282376464804295311802328782509819455815301756717361332069811250996181881593041690351598888519345807273866738589422879228499892086805825749279610484198444363463244968487560233624827041978623209002160990235304369941849146314093431738143640546253152096183690888707016768396424378140592714563549061303107208510383750510115747704171898610687396965521267154688957035035";
}