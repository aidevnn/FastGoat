using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Integers;

public enum PadicDisplayForm
{
    Serie,
    UnsignedInteger,
    SignedInteger,
    Number
}

public readonly struct Padic : IEquatable<Padic>, IComparable<Padic>
{
    public static PadicDisplayForm DisplayForm { get; } = PadicDisplayForm.Number;
    private Dictionary<int, int> Digits { get; }
    public int P { get; }
    public int O { get; }
    public int Emin => Digits.Count == 0 ? 0 : Int32.Min(0, Digits.Keys.Min());
    public int Emax => Digits.Count == 0 ? 0 : Int32.Max(0, Digits.Keys.Max());

    public Padic()
    {
        Digits = new();
        P = 2;
        O = 1;
    }

    public Padic(int p, int o)
    {
        Digits = new();
        P = p;
        O = o;
    }

    private Padic((int p, int o) e, int capacity)
    {
        Digits = new(capacity);
        P = e.p;
        O = e.o;
    }

    public Padic(int p, int o, BigInteger v)
    {
        var p0 = Convert(p, v);
        Digits = p0.Digits;
        P = p;
        O = o;
    }

    public Padic(Padic c)
    {
        Digits = new(c.Digits);
        P = c.P;
        O = c.O;
    }

    public Padic Resize(int o)
    {
        var a = new Padic(P, o);
        for (int i = Emin; i <= o; i++)
            a[i] = this[i];
        return a;
    }

    public int Valuation => -FirstCoeff().r;

    public (int r, int Coeff) FirstCoeff()
    {
        for (int r = Emin; r < O; ++r)
        {
            var v = this[r];
            if (v != 0)
                return (r, v);
        }

        return (0, 0);
    }

    public Padic Normalized => Shift(Valuation);

    public (int E, int H, Padic Sci) ScientificForm()
    {
        var e = Digits.Count == 0 ? 0 : Digits.Keys.Min();
        var a = new Padic(P, O - e);
        for (int i = e; i < O; i++)
            a[i - e] = this[i];

        var h = Digits.Count == 0 ? 0 : a.Digits.Keys.Max();
        return (e, h, a);
    }

    public bool IsInteger() => Emin == 0;

    public Padic IntegerPart
    {
        get
        {
            var a = new Padic((P, O), O);
            for (int i = 0; i < O; i++)
                a[i] = this[i];

            return a;
        }
    }

    public Padic FloatingPart
    {
        get
        {
            var a = new Padic((P, O), O);
            for (int i = Emin; i < 0; i++)
                a[i] = this[i];

            return a;
        }
    }

    public Padic Shift(int k)
    {
        var a = new Padic(P, O);
        for (int i = Emin; i < O; i++)
            a[i - k] = this[i];

        return a;
    }

    public BigInteger ToBigInteger()
    {
        var (e, _, s) = ScientificForm();
        var acc = BigInteger.Zero;
        var pow = BigInteger.Pow(P, e);
        for (int i = 0; i < s.O; i++)
        {
            acc += s[i] * pow;
            pow *= P;
        }

        return acc;
    }

    public bool IsZero() => Digits.Count == 0;

    public Padic Zero => new(P, O);

    public Padic One => new(P, O, 1);

    public Padic Add(Padic other)
    {
        if (P != other.P || O != other.O)
            throw new ArgumentException();

        var e = Int32.Min(Emin, other.Emin);
        var add = new Padic((P, O), O - e);
        var carry = 0;
        for (int i = e; i < O; i++)
        {
            var s = this[i] + other[i] + carry;
            if (s < P)
            {
                add[i] = s;
                carry = 0;
            }
            else
            {
                add[i] = s - P;
                carry = 1;
            }
        }

        return add;
    }

    public Padic Sub(Padic other)
    {
        if (P != other.P || O != other.O)
            throw new ArgumentException();

        var e = Int32.Min(Emin, other.Emin);
        var sub = new Padic((P, O), O - e);
        var carry = 0;
        for (int i = e; i < O; i++)
        {
            var s = this[i] - (other[i] + carry);
            if (s >= 0)
            {
                sub[i] = s;
                carry = 0;
            }
            else
            {
                sub[i] = s + P;
                carry = 1;
            }
        }

        return sub;
    }

    public Padic Opp() => Zero.Sub(this);

    public Padic Mul(Padic other)
    {
        if (P != other.P || O != other.O)
            throw new ArgumentException();

        var e = Emin + other.Emin;
        var mul = new Padic((P, O), O - e);
        var carry = 0;
        for (int k = e; k < O; k++)
        {
            var sum = carry;
            for (int i = Emin; i < O; i++)
                sum += this[i] * other[k - i];

            var r = mul[k] = sum % P;
            carry = (sum - r) / P;
        }

        return mul;
    }

    public Padic Mul(int k)
    {
        var e = new Padic(P, O, k);
        return Mul(e);
    }

    public (Padic quo, Padic rem) Div(Padic other)
    {
        if (P != other.P || O != other.O)
            throw new ArgumentException();

        var a0 = CompleteDivision(other);
        var quo = a0.IntegerPart;
        var rem = Sub(quo.Mul(other));
        return (quo, rem);
    }

    public Padic CompleteDivision(Padic other)
    {
        if (P != other.P || O != other.O)
            throw new ArgumentException();

        var (v1, h1, n1) = other.ScientificForm();
        var o1 = O + Int32.Abs(v1) + 2;
        n1 = n1.Resize(o1);
        var a0 = new Padic(this).Resize(o1);
        var quo = new Padic(P, n1.O);

        var e1 = IntExt.InvModP(n1[0], P);
        while (!a0.IsZero())
        {
            var (v0, n0) = a0.FirstCoeff();
            var c0 = (n0 * e1) % P;
            quo[v0] = c0;
            a0 = a0.SubMul(n1, c0, v0);
            // a0 = a0.Sub(n1.Mul(c0).Shift(v0));
        }

        return quo.Shift(v1).Resize(O);
    }

    private Padic SubMul(Padic other, int k, int s)
    {
        var submul = new Padic(P, O);
        var e = Emin;
        var carry = 0;
        for (int i = e; i < O; i++)
        {
            var sum = this[i] - (other[i - s] * k + carry);
            var r0 = sum % P;
            var r = submul[i] = r0 < 0 ? P + r0 : r0;
            var c0 = (sum - r) / P;
            carry = c0 < 0 ? -c0 : c0;
        }

        return submul;
    }

    public Padic Inv() => One.CompleteDivision(this);

    public bool Equals(Padic other)
    {
        if (P != other.P || O != other.O || Emin != other.Emin)
            return false;

        for (int i = Emin; i < O; i++)
        {
            if (this[i] != other[i])
                return false;
        }

        return true;
    }

    public int CompareTo(Padic other)
    {
        if (P != other.P || O != other.O)
            throw new ArgumentException();

        var e = Int32.Min(Emin, other.Emin);
        for (int i = O - 1; i >= e; i--)
        {
            var comp = this[i].CompareTo(other[i]);
            if (comp != 0)
                return comp;
        }

        return 0;
    }

    public override int GetHashCode()
    {
        var acc = P;
        for (int i = Emin; i < O; i++)
            acc = (acc, this[i]).GetHashCode();

        return acc;
    }

    public int this[int i]
    {
        get
        {
            if (Digits.ContainsKey(i))
                return Digits[i];

            return 0;
        }
        private set
        {
            if (i > O + 2)
                return;

            if (value == 0)
                Digits.Remove(i);
            else
                Digits[i] = value;
        }
    }

    public override string ToString()
    {
        var pstr = Ring.Xi('p', O).ToString().Replace("p", $"{P}");
        var s = "";

        if (IsZero())
            return $"[0({pstr})]";

        for (int i = Emin; i <= Int32.Min(O - 1, Emax); i++)
        {
            var sep = P < 10 ? "" : "'";
            sep = i == 1 ? "." : sep;
            sep = i == Emin ? "" : sep;
            s = $"{this[i]}" + sep + s;
        }

        s = s.Reverse().Glue();
        return $"[{s}({pstr})]";
    }

    public static Padic Convert(int p, BigInteger v)
    {
        var Digits = new Dictionary<int, int>() { [0] = 0 };
        var v0 = v < 0 ? -v : v;
        var i = 0;
        while (v0 != 0)
        {
            var r = v0 % p;
            Digits[i] = r < 0 ? (int)(r + p) : (int)r;
            v0 = (v0 - r) / p;
            ++i;
        }

        var o = Digits.Keys.Max();
        var a = new Padic(p, o + 1);
        for (int j = 0; j <= o; j++)
            a[j] = Digits[j];

        return v < 0 ? a.Opp() : a;
    }

    public static Padic Convert(int p, int o, Rational r)
    {
        var num0 = r.Num >= 0 ? r.Num : -r.Num;
        var denom = Convert(p, r.Denom);
        var o1 = Int32.Max(denom.O, o);
        var num = new Padic(p, o1, num0);
        denom = denom.Resize(o1);
        Console.WriteLine($"num={r.Num.Sign}*{num} denom={denom}");
        if (r.Num >= 0)
            return num.CompleteDivision(denom).Resize(o);

        return num.CompleteDivision(denom).Opp().Resize(o);
    }

    public static Padic Convert(int p, int o, (BigInteger num, BigInteger denom) e) =>
        Convert(p, o, new Rational(e.num, e.denom));
}