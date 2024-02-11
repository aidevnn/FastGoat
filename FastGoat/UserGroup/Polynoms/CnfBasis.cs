using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Polynoms;

public struct CnfBasis
{
    public int Order { get; }
    public int[] Factors { get; }
    public KPoly<Rational> CyclotomicPolynomial { get; }
    public SortedList<int, KMatrix<Rational>> Basis { get; }
    public Dictionary<int, int[]> RelativeIndexes { get; }
    public ConcreteGroup<Ep<ZnInt>> Gr { get; }

    public Dictionary<(int, bool), Cnf> IdxRoots { get; }

    public CnfBasis(int ord)
    {
        Order = ord;
        CyclotomicPolynomial = FG.CyclotomicPolynomial(ord);
        var deg = CyclotomicPolynomial.Degree;
        var x = FG.EPoly(CyclotomicPolynomial, Cnf.RootsOfUnit);
        Factors = IntExt.PrimesDec(Order).Select(kv => kv.Key.Pow(kv.Value)).Order().ToArray();
        IdxRoots = ord.Range().Select(k => (Cnf.Nth(ord).Pow(k), k))
            .SelectMany(e => new[] { (e.Item1, (e.k, true)), (-e.Item1, (e.k, false)) })
            .ToDictionary(e => e.Item2, e => e.Item1);
        var gr = Gr = FG.Abelian(Factors);
        var allSubgroups = Group.AllSubGroups(Gr).Keys.Where(sg => sg.Count() != 1).ToHashSet(new IsomorphEquality<Ep<ZnInt>>());

        RelativeIndexes = new(Factors.Length);
        Basis = new(Factors.Length);
        foreach (var sg in allSubgroups)
        {
            var key = sg.Count();
            var elts = sg.OrderByDescending(e => gr.ElementsOrders[e]).ToArray();

            KMatrix<Rational> bs = new(Rational.KZero(), deg, deg); // useless initialization
            var idxs = new List<int>();
            var start = true;
            foreach (var elt in elts)
            {
                var idxRel = elt.Ei.Select(z => (z.K * key / z.P) % key).Sum() % key;
                var idxAbs = elt.Ei.Select(z => (z.K * ord / z.P) % ord).Sum() % ord;
                if (start)
                {
                    bs = x.Pow(idxAbs).Poly.ToVMatrix(deg);
                    idxs.Add(idxRel);
                    start = false;
                    continue;
                }

                var mat = KMatrix<Rational>.MergeSameRows(bs, x.Pow(idxAbs).Poly.ToVMatrix(deg));
                var (dim, ker) = mat.NullSpace();
                if (dim == 0)
                {
                    bs = mat;
                    idxs.Add(idxRel);
                }
            }

            Basis[key] = bs;
            RelativeIndexes[key] = idxs.ToArray();
        }
    }

    public (Cnf, KPoly<Rational>) Simplify(Cnf c)
    {
        if (c.N == 1 || c.N == 4)
        {
            var p = c.E.Poly.Div(c.E.F).rem;
            return (new(c.N, new(CyclotomicPolynomial, p)), p);
        }

        if (c.E.Degree == 0)
        {
            var p = c.E.Poly;
            return (new(1, new(p.X - 1, p)), p);
        }

        if (c.N != Order)
            throw new();

        var matC = c.E.Poly.ToVMatrix(CyclotomicPolynomial.Degree);
        foreach (var (ord, matBs) in Basis)
        {
            var mat = KMatrix<Rational>.MergeSameRows(matBs, matC);
            var (dim, ker) = mat.NullSpace();

            if (dim == 1)
            {
                var idx = RelativeIndexes[ord];
                var x = FG.CyclotomicEPoly(ord);
                var arr = (ord + 1).Range().Select(_ => Rational.KZero()).ToArray();
                foreach (var (k, v) in idx.Select((k, i) => (k, -ker[i, 0])))
                    arr[k] = v;

                var p0 = new KPoly<Rational>('x', Rational.KZero(), arr.TrimSeq().ToArray());
                var c0 = new Cnf(ord, new(x.F, p0.Div(x.F).rem));
                if (arr.Count(r => !r.IsZero()) > 1)
                {
                    var (re, im) = (c.Re, c.Im);
                    var mod2 = re * re + im * im;

                    if (mod2.Equals(Cnf.CnfOne) && !re.IsZero() && !im.IsZero())
                    {
                        var (k, a) = IdxRoots.First(e => (e.Value.E - c.E).IsZero()).Key;
                        var gcd = IntExt.Gcd(k, Order);
                        var (k0, n0) = (k / gcd, Order / gcd);
                        var c1 = Cnf.Nth(n0).Pow(k0);
                        var p1 = FG.QPoly().Pow(k0);
                        if (!a)
                        {
                            c1 = -c1;
                            p1 = -p1;
                        }

                        return (c1, p1);
                    }
                }

                return (c0, p0);
            }
        }

        throw new();
    }
}