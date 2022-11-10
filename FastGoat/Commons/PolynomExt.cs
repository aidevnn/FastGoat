namespace FastGoat.Commons;

public static class PolynomExt
{
    static PolynomExt()
    {
        All = new();
        foreach (var str in Table)
        {
            var split = str.Split(':');
            var pm = split[0].Split(',').Select(int.Parse).ToArray();
            var coefs = split[1].Split(',').Select(int.Parse).ToArray();

            var p = pm[0];
            var m = pm[1];
            if (!All.ContainsKey(p))
                All[p] = new();

            All[p][m] = coefs;
        }
    }

    private static Dictionary<int, Dictionary<int, int[]>> All;

    public static ((int p, int m) pm, int[] coefs) Get(int q)
    {
        var pm = IntExt.PrimesDecomposition(q).ToArray();
        if (pm.Distinct().Count() != 1)
            throw new();

        var p = pm[0];
        var m = pm.Length;
        return ((p, m), All[p][m]);
    }

    public static string PolyStr(char x, string prev, int n, int v)
    {
        if (v == 0) return prev;

        int sgn = v < 0 ? -1 : 1;
        int va = v < 0 ? -v : v;
        var xs = n == 0 ? "" : n == 1 ? $"{x}" : $"{x}^{n}";
        var vs = n == 0 ? $"{va}" : va == 1 ? $"{xs}" : $"{va}{xs}";
        if (prev == "")
            return sgn == -1 ? $"-{vs}" : $"{vs}";

        return sgn == -1 ? $"{prev} - {vs}" : $"{prev} + {vs}";
    }

    public static int[] TrimPoly(int[] poly)
    {
        var stack = new Stack<int>(poly);
        while (stack.Count != 0 && stack.Peek() == 0)
            stack.Pop();

        if (stack.Count == 0)
            return new[] { 0 };

        return stack.Reverse().ToArray();
    }

    static void AddInPlace(int p, int start, int[] acc, int[] m, int factor = 1)
    {
        for (int i = 0; i < m.Length; i++)
            acc[start + i] = IntExt.AmodP(acc[start + i] + m[i] * factor, p);
    }

    public static int[] AddPoly(int p, int[] a, int[] b)
    {
        if (a.Length >= b.Length)
        {
            var acc = a.ToArray();
            AddInPlace(p, 0, acc, b);
            return TrimPoly(acc);
        }
        else
        {
            var acc = b.ToArray();
            AddInPlace(p, 0, acc, a);
            return TrimPoly(acc);
        }
    }

    public static int[] MulKPoly(int p, int k, int[] b)
    {
        return TrimPoly(b.Select(e => IntExt.AmodP(e * k, p)).ToArray());
    }

    public static int[] MulPoly(int p, int[] a, int[] b)
    {
        var acc = new int[a.Length + b.Length];
        for (int i = 0; i < a.Length; i++)
            AddInPlace(p, i, acc, b, a[i]);

        return TrimPoly(acc);
    }

    public static (int[] quo, int[] rem) DivPoly(int p, int[] a, int[] b)
    {
        if (a.Length >= b.Length)
        {
            var bn = b.Last();
            var bi = IntExt.InvModP(bn, p);
            var quo = new int[a.Length - b.Length + 1];
            var rem = a.ToArray();
            for (int i = 0, k = a.Length - 1; i < quo.Length; ++i, --k)
            {
                var f = IntExt.AmodP(-rem[k] * bi, p);
                quo[quo.Length - i - 1] = IntExt.AmodP(-f, p);
                AddInPlace(p, quo.Length - i - 1, rem, b, f);
            }

            return (quo, TrimPoly(rem));
        }
        else
        {
            return (new[] { 0 }, a);
        }
    }

    // http://www.math.rwth-aachen.de/~Frank.Luebeck/data/ConwayPol/index.html
    static readonly string[] Table = new[]
    {
        "2,1:1,1",
        "2,2:1,1,1",
        "2,3:1,1,0,1",
        "2,4:1,1,0,0,1",
        "2,5:1,0,1,0,0,1",
        "2,6:1,1,0,1,1,0,1",
        "2,7:1,1,0,0,0,0,0,1",
        "2,8:1,0,1,1,1,0,0,0,1",
        "2,9:1,0,0,0,1,0,0,0,0,1",
        "3,1:1,1",
        "3,2:2,2,1",
        "3,3:1,2,0,1",
        "3,4:2,0,0,2,1",
        "3,5:1,2,0,0,0,1",
        "3,6:2,2,1,0,2,0,1",
        "3,7:1,0,2,0,0,0,0,1",
        "3,8:2,2,2,0,1,2,0,0,1",
        "3,9:1,1,2,2,0,0,0,0,0,1",
        "5,1:3,1",
        "5,2:2,4,1",
        "5,3:3,3,0,1",
        "5,4:2,4,4,0,1",
        "5,5:3,4,0,0,0,1",
        "5,6:2,0,1,4,1,0,1",
        "5,7:3,3,0,0,0,0,0,1",
        "5,8:2,4,3,0,1,0,0,0,1",
        "5,9:3,1,0,2,0,0,0,0,0,1",
        "7,1:4,1",
        "7,2:3,6,1",
        "7,3:4,0,6,1",
        "7,4:3,4,5,0,1",
        "7,5:4,1,0,0,0,1",
        "7,6:3,6,4,5,1,0,1",
        "7,7:4,6,0,0,0,0,0,1",
        "7,8:3,2,6,4,0,0,0,0,1",
        "7,9:4,6,0,1,6,0,0,0,0,1",
        "11,1:9,1",
        "11,2:2,7,1",
        "11,3:9,2,0,1",
        "11,4:2,10,8,0,1",
        "11,5:9,0,10,0,0,1",
        "11,6:2,7,6,4,3,0,1",
        "11,7:9,4,0,0,0,0,0,1",
        "11,8:2,7,1,7,7,0,0,0,1",
        "11,9:9,8,9,0,0,0,0,0,0,1",
        "13,1:11,1",
        "13,2:2,12,1",
        "13,3:11,2,0,1",
        "13,4:2,12,3,0,1",
        "13,5:11,4,0,0,0,1",
        "13,6:2,11,11,10,0,0,1",
        "13,7:11,3,0,0,0,0,0,1",
        "13,8:2,3,2,12,8,0,0,0,1",
        "13,9:11,12,12,8,12,0,0,0,0,1",
        "17,1:14,1",
        "17,2:3,16,1",
        "17,3:14,1,0,1",
        "17,4:3,10,7,0,1",
        "17,5:14,1,0,0,0,1",
        "17,6:3,3,10,0,2,0,1",
        "17,7:14,12,0,0,0,0,0,1",
        "17,8:3,6,0,12,11,0,0,0,1",
        "17,9:14,8,7,0,0,0,0,0,0,1",
        "19,1:17,1",
        "19,2:2,18,1",
        "19,3:17,4,0,1",
        "19,4:2,11,2,0,1",
        "19,5:17,5,0,0,0,1",
        "19,6:2,6,17,17,0,0,1",
        "19,7:17,6,0,0,0,0,0,1",
        "19,8:2,3,10,12,1,0,0,0,1",
        "19,9:17,16,14,11,0,0,0,0,0,1",
        "23,1:18,1",
        "23,2:5,21,1",
        "23,3:18,2,0,1",
        "23,4:5,19,3,0,1",
        "23,5:18,3,0,0,0,1",
        "23,6:5,1,9,9,1,0,1",
        "23,7:18,21,0,0,0,0,0,1",
        "23,8:5,3,5,20,3,0,0,0,1",
        "23,9:18,9,8,3,0,0,0,0,0,1"
    };
}