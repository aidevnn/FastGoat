using System.Numerics;
using System.Reflection;
using System.Text;
using Craft;
using Craft.Craft;
using Examples;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.Tools;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

Perm.Style = DisplayPerm.CyclesComma;

ConcreteGroup<Perm> ProductPermGroup(params ConcreteGroup<Perm>[] Gs)
{
    var HnBase = Product.Gp(Gs.Cast<IGroup<Perm>>().ToArray());
    var HnGens = HnBase.GetGenerators().Select(e => FG.ConcatPerm(e.Ei)).ToArray();
    return Group.Generate(HnBase.Name, HnGens[0].Sn, HnGens);
}

bool TypeMatch(ConcreteGroup<Automorphism<Perm>> autG)
{
    var G = autG.Neutral().Domain;
    return autG.GetGenerators().All(aut => G.GetGenerators().All(e => Perm.TypeEquals(e, aut[e])));
}

string GapExport(Perm[] gens)
{
    var old = Perm.Style;
    Perm.Style = DisplayPerm.Gap;
    var export = $"Group([{gens.Glue(", ")}]);";
    Perm.Style = old;
    return export;
}

ConcreteGroup<Automorphism<T>> AutomorphismGroup<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
{
    var autbs = new AutomorphismGroup<T>(g);
    var allAut = GroupCraft.AllMorphismsWithPruning(g, g, Group.MorphismType.Isomorphism)
        .Select(aut => new Automorphism<T>(autbs, aut.HomMap))
        .OrderByDescending(aut => Group.GenerateElements(autbs, aut).Count)
        .ToArray();
    var autG = Group.Generate($"Aut[{g.Name}]", autbs, allAut);
    return autG;
}

ConcreteGroup<Perm> RegPermAutGroup<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
{
    var autbs = new AutomorphismGroup<T>(G);
    var sn = new Sn(G.Count());
    var mapEltIdx = G.OrderBy(e => G.ElementsOrders[e]).Index().ToDictionary(e => e.Item, e => e.Index);
    var allAut = GroupCraft.AllMorphismsWithPruning(G, G, Group.MorphismType.Isomorphism)
        .Select(aut => new Automorphism<T>(autbs, aut.HomMap))
        .Select(e => e.AutMap.ToDictionary(f => mapEltIdx[f.Key], f => mapEltIdx[f.Value]))
        .Select(e => sn.CreateElementTable(e.OrderBy(f => f.Key).Select(f => f.Value).ToArray()))
        .ToArray();
    var autG = Group.Generate($"Aut[{G.Name}]", sn, allAut);
    return autG;
}

ConcreteGroup<Perm> RegPermZentrumAutGroup<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
{
    var autG = RegPermAutGroup(G);
    var sn = autG.Neutral().Sn;
    var actSn = Group.ByConjugate(sn);
    var zentrumAut = autG.Where(e => autG.GetGenerators().All(aut => actSn.IsInvariant(e, aut))).ToArray();
    var zAutG = Group.Generate($"Z[Aut[{G.Name}]]", sn, zentrumAut);
    var zType = Group.AbelianDecompositions(zAutG).abType.Select(e => e.g).ToArray();
    return Group.Generate($"Z[Aut[{G.Name}]]", sn, zType);
}

ConcreteGroup<Perm> RegPermInnAutGroup<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
{
    var sn = new Sn(G.Count());
    var mapEltIdx = G.OrderBy(e => G.ElementsOrders[e]).Index().ToDictionary(e => e.Item, e => e.Index);
    var act = Group.ByConjugate(G);
    var autbs = Group.AutBase(G);
    var gensInnAut = G.GetGenerators().Select(e => G.GetGenerators().ToDictionary(f => f, f => act(e, f)))
        .Select(pmap => new Automorphism<T>(autbs, Group.IsomorphismMap(G, G, pmap)))
        .Select(e => e.AutMap.ToDictionary(f => mapEltIdx[f.Key], f => mapEltIdx[f.Value]))
        .Select(e => sn.CreateElementTable(e.OrderBy(f => f.Key).Select(f => f.Value).ToArray()))
        .ToArray();

    return Group.Generate($"InnAut[{G.Name}]", sn, gensInnAut);
}

(ConcreteGroup<Perm> G, ConcreteGroup<Perm> AutG, Dictionary<T, Perm> MapReg)
    RegPermGroupAndAutGroup<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
{
    var sn = new Sn(G.Count());
    var mapEltIdx = G.Index().OrderBy(e => G.ElementsOrders[e.Item])
        .ToDictionary(e => e.Item, e => e.Index);
    var mapIdxElt = mapEltIdx.ToDictionary(e => e.Value, e => e.Key);
    GroupAction<T, Perm> act = (e, p) => sn.CreateElementTable(p.Table.Select(i => mapEltIdx[G.Op(e, mapIdxElt[i])]).ToArray());
    var mapReg = G.ToDictionary(e => e, e => act(e, sn.Neutral()));
    var gensReg = G.GetGenerators().Select(e => mapReg[e]).ToArray();
    var Greg = Group.Generate(G.Name, sn, gensReg);

    var autbs = new AutomorphismGroup<T>(G);
    var allAut = GroupCraft.AllMorphismsWithPruning(G, G, Group.MorphismType.Isomorphism)
        .Select(aut => new Automorphism<T>(autbs, aut.HomMap))
        .Select(e => e.AutMap.ToDictionary(f => mapEltIdx[f.Key], f => mapEltIdx[f.Value]))
        .Select(e => sn.CreateElementTable(e.OrderBy(f => f.Key).Select(f => f.Value).ToArray()))
        .OrderByDescending(e => e.Order)
        .ToArray();
    var autG = Group.Generate($"Aut[{G.Name}]", sn, allAut);
    return (Greg, autG, mapReg);
}

void DisplaySDPAb<T>(AllSubgroups<T> subs) where T : struct, IElt<T>
{
    subs.DecomposeProducts(subs.ProperNonTrivialNormalSubgroups())
        .Where(e => e.lhs.GroupType == GroupType.AbelianGroup && e.rhs.GroupType == GroupType.AbelianGroup)
        .Select(e => (lhs: Group.AbelianGroupType(e.lhs.Representative),
            rhs: Group.AbelianGroupType(e.rhs.Representative),
            symb: e.isDirectProduct ? "x" : "x:"))
        .DistinctBy(e => e.lhs.ToAbString())
        .Println(e => $"{e.lhs.ToAbString().WithParenthesis()} {e.symb} {e.rhs.ToAbString().WithParenthesis()}",
            subs.Parent.Name);
    Console.WriteLine();
}

void DisplaySDPAbLhs<T>(AllSubgroups<T> subs, int[] zType) where T : struct, IElt<T>
{
    var seq = subs.Where(e => e.IsNormal && e.GroupType == GroupType.AbelianGroup)
        .Select(e => (lhs: e, lhsType: Group.AbelianGroupType(e.Representative)))
        .Where(e => AbelianExt.IsSubgroup(e.lhsType, zType))
        .Select(e => (e.lhs,
            e.lhsType,
            comps: subs.Complements(e.lhs)
                .DistinctBy(r => r.fact.Representative.ElementsOrdersList().Glue(",")).ToArray()))
        .SelectMany(e => e.comps.Select(rhs => (e.lhs, e.lhsType, rhs: rhs.fact, rhs.isDirectprod)))
        .DistinctBy(e => $"{e.lhsType.Glue(",")} {e.rhs.Representative.ElementsOrdersList().Glue(",")}")
        .ToArray();
    Console.WriteLine(subs.Parent.Name);
    foreach (var (lhs, lhsType, rhs, isDirectprod) in seq)
    {
        var symb = isDirectprod ? "x" : "x:";
        var rhsRes = subs.Restriction(rhs.Representative);
        rhsRes.Naming();
        var rhsName = rhsRes.Parent.Name;
        var diff = lhsType.SequenceEqual(zType) ? "zType" : "superset";
        Console.WriteLine(
            $"    {($"{lhsType.ToAbString().WithParenthesis()} {symb} {rhsName.WithParenthesis()}"),-60} diff:{diff}");
    }

    if (seq.Length == 0)
        Console.WriteLine("    empty");

    Console.WriteLine();
}

Perm[] CyclesSplit(int m)
{
    var dec = IntExt.PrimesDec(m).Select(e => e.Key.Pow(e.Value)).Order().ToArray();
    return dec.Select(e => new Sn(e).Cycle(e.Range(1))).ToArray();
}

Perm rUmToPerm(int m, int r)
{
    var cycles = CyclesSplit(m).Select(c =>
    {
        var N = c.Sn.N;
        if (IntExt.Gcd(N, r) != 1)
            return c;

        var ri = new ZnInt(N, r).Inv();
        return c.Sn.CreateElementTable(N.Range().Select(i => (i * ri).K).ToArray());
    }).ToArray();

    return FG.ConcatPerm(cycles.ToArray());
}

Perm[] rUmToPermAll(int m, int r)
{
    var a = FG.Cycles(m);
    var b = a ^ r;
    var sn = a.Sn;
    var aOrbits = a.Orbits.OrderBy(e => e.Length).ToArray();
    var bOrbits = b.Orbits.OrderBy(e => e.Length).ToArray();
    return bOrbits.Select(c => UGCraft.CycleWalk(c)).MultiLoop().Select(l => l.ToArray())
        .Select(cycle =>
        {
            var autTable = new int[sn.N];
            foreach (var (bOrb, aOrb) in cycle.Zip(aOrbits))
                for (int i = 0; i < bOrb.Length; i++)
                    autTable[bOrb[i]] = aOrb[i];

            return sn.CreateElementTable(autTable);
        }).ToArray();
}

(string name, Sn n, Perm a, Perm b) MetaCyclicGens(int m, int n, int r)
{
    var a = FG.Cycles(m);
    var b0 = rUmToPerm(m, r);
    var b = FG.Cycles(n);
    var name = FG.MetaCyclicName(m, n, r);
    return (name, new(a.Sn.N + b.Sn.N), FG.PaddingRight(a, b.Sn.N), FG.ConcatPerm(b0, b));
}

void OrderedMetaCyclic(int maxOrder)
{
    var seqm = (maxOrder / 2 - 2).Range(3).ToDictionary(m => m, m => IntExt.Coprimes(m).Where(u => u > 1).ToXSet());
    var seqn = (maxOrder / 3 - 1).Range(2).ToDictionary(n => n, n => IntExt.Coprimes(n).ToXSet());
    var all = seqm.Grid2D(seqn)
        .SelectMany(e => e.t1.Value.Where(u => IntExt.PowMod(u, e.t2.Key, e.t1.Key) == 1)
            .OrderDescending()
            .Select(r => (e.t1, e.t2, r)))
        .DistinctBy(e => (e.t1.Key, e.t2.Key, e.t2.Value.Select(i => IntExt.PowMod(e.r, i, e.t1.Key)).ToXSet()))
        .Select(e => (m: e.t1.Key, n: e.t2.Key, e.r, mCoprimes: e.t1.Value, nCoprimes: e.t2.Value))
        .Where(e => e.m * e.n <= maxOrder)
        .GroupBy(e => (e.m, e.r))
        .ToDictionary(gr => gr.Key, gr => gr.ToArray());

    foreach (var ((m, r), gr) in all.OrderBy(e => e.Key))
        gr.Println(e => $"   n:{e.n,-4} {FG.MetaCyclicName(m, e.n, e.r)}", $"m:{m,-4} r:{r}");

    // foreach (var gr in all.Values.SelectMany(e => e)
    //              .GroupBy(e => (e.m, e.n))
    //              .OrderBy(e => (e.Key.m, e.Key.n)))
    // {
    //     gr.Order().Println(e => FG.MetaCyclicName(e.m, e.n, e.r), $"m:{gr.Key.m} n:{gr.Key.n}");
    // }
}

Dictionary<int, (int m, int n, int r)[]> MetaCyclicsM(int m, int maxOrder, bool details = false)
{
    var maxN = maxOrder / m;
    var mCoprimes = IntExt.Coprimes(m).Where(u => u > 1).ToXSet();
    var seqn = (maxN - 1).Range(2).ToDictionary(n => n, n => IntExt.Coprimes(n).ToXSet());
    var all = seqn.SelectMany(e => mCoprimes.Where(u => IntExt.PowMod(u, e.Key, m) == 1)
            .OrderDescending()
            .Select(r => (m, e, r)))
        .DistinctBy(e => (e.e.Key, e.e.Value.Select(i => IntExt.PowMod(e.r, i, m)).ToXSet()))
        .Select(e => (n: e.e.Key, e.r, mCoprimes, nCoprimes: e.e.Value))
        .Where(e => m * e.n <= maxOrder)
        .GroupBy(e => e.r)
        .OrderBy(e => e.Key)
        .ToDictionary(gr => gr.Key, gr => gr.ToArray());

    if (details)
    {
        foreach (var (r, gr) in all)
            gr.Println(e => $"   n:{e.n,-4} {FG.MetaCyclicName(m, e.n, e.r)}", $"m:{m,-4} r:{r}");
    }

    return all.ToDictionary(e => e.Key, e => e.Value.Select(f => (m, f.n, f.r)).ToArray());
}

(ConcreteGroup<Perm> mtp, ConcreteGroup<Perm> autMtp) AutMetaCyclicPNRSDP(int p, int n, int r)
{
    if (IntExt.PrimesDec(p).Count != 1 || IntExt.Gcd(p, n) != 1)
        throw new();

    if (r == 1)
    {
        var (ab, autAb) = (FG.AbelianPerm(p * n), FG.Dihedral(p * n));
        DisplayGroup.HeadOrders(ab);
        autAb.Name = $"Aut[{ab}] = {autAb}";
        DisplayGroup.HeadOrders(autAb);
        return (ab, autAb);
    }

    var mtp = FG.MetaCyclicPg(p, n, r);
    DisplayGroup.HeadOrdersGenerators(mtp);

    var q = IntExt.Phi(p);
    var rq = IntExt.Coprimes(p).OrderDescending().First(rq => rUmToPerm(p, rq).Order == q);
    var mtpq = FG.MetaCyclicPg(p, q, rq);

    var unType = Group.AbelianGroupType(new Un(n));
    var pr = rUmToPerm(p, r).Order;
    var ubType = Group.AbelianGroupType(new Un(pr));
    var quos = AbelianExt.QuotientsCan(unType, ubType);
    if (quos.Length > 1)
    {
        quos.Println(l => l.ToAbString(), $"U{n}={unType.ToAbString()} U{pr}:{ubType.ToAbString()}");
        Console.WriteLine("Warning"); // TODO: quotient criterion
    }

    var lhsType = quos.MinBy(l => l.Length)!;
    var lhsAb = FG.AbelianPerm(lhsType);
    var autMtpg = ProductPermGroup(mtpq, lhsAb);
    autMtpg.Name = $"Aut[{mtp}]";
    Console.WriteLine($"{autMtpg} = {lhsAb.Name} x {mtpq.Name}");
    Console.WriteLine();
    DisplayGroup.HeadOrdersGenerators(autMtpg);

    if (autMtpg.Order <= 1024)
    {
        var autMt = RegPermGroupAndAutGroup(mtp).AutG;
        DisplayGroup.HeadOrders(autMt);
        Console.WriteLine();
        if (autMtpg.ElementsOrdersList().SequenceEqual(autMt.ElementsOrdersList()))
            return (mtp, autMtpg);

        throw new();
    }

    Console.WriteLine($"mtp:={GapExport(mtp.GetGenerators().ToArray())}");
    Console.WriteLine($"autMt:={GapExport(autMtpg.GetGenerators().ToArray())}");
    Console.WriteLine("StructureDescription(autMt);");
    Console.WriteLine();

    return (mtp, autMtpg);
}

void AutMetacyclicP(int p, int maxOrder)
{
    if (!IntExt.IsPrime(p))
        throw new();
    var mtpCoefs = MetaCyclicsM(p, maxOrder, details: true);
    var (counter, found, check) = (0, 0, 0);
    foreach (var (_, coefs) in mtpCoefs)
    {
        foreach (var (_, n, r) in coefs)
        {
            ++counter;
            // Aut[M(p,n)r] = M x Fpq where p is prime, q=p-1 and M is an abelian group
            var (mtp, autMtp) = AutMetaCyclicPNRSDP(p, n, r);
            if (autMtp.Order != 1)
                ++found;
        }
    }

    Console.WriteLine($"Total:{counter} Found:{found}");
    Console.WriteLine();
}

void AllAutMetacyclicP(int maxP, int maxOrder)
{
    var counter = 0;
    foreach (var p in IntExt.Primes10000.Where(p => p <= maxP))
    {
        var mtpCoefs = MetaCyclicsM(p, maxOrder, details: true);
        foreach (var (_, coefs) in mtpCoefs)
        {
            foreach (var (_, n, r) in coefs)
            {
                ++counter;
                // Aut[M(p,n)r] = M x Fpq where p is prime, q=p-1 and M is an abelian group
                AutMetaCyclicPNRSDP(p, n, r);
            }
        }
    }

    Console.WriteLine($"Total:{counter}");
    Console.WriteLine();
}

(string name, int[] abType, int m, int n, int r) RewriteMetaCyclic(int m, int n, int r)
{
    var decM = IntExt.PrimesDec(m).Select(e => e.Key.Pow(e.Value)).ToArray();
    var mAbType = decM.Where(e => r % e == 1).ToArray();
    var m2 = m / mAbType.Aggregate(1, (acc, e) => e * acc);
    var r2 = r % m2;

    var decN = IntExt.PrimesDec(n).Select(e => e.Key.Pow(e.Value)).ToArray();
    var nAbType = decN.Where(e => IntExt.PowMod(r2, e, m2) != 1).ToArray();
    var n2 = n / nAbType.Aggregate(1, (acc, e) => e * acc);
    if (n2 == 1)
        (n2, nAbType) = (n, []);

    var abType = AbelianExt.AbType(mAbType.Concat(nAbType).ToArray()).can;
    var name = FG.MetaCyclicName(m2, n2, r2);
    if (abType[0] != 1)
        name = $"{abType.ToAbString()} x {name}";
    return (name, abType, m2, n2, r2);
}

void TestDecomposeMetaCyclic(int maxM, int maxOrder)
{
    var (total, found, notFound) = (0, 0, 0);
    foreach (var m in (maxM - 2).SeqLazy(3))
    {
        var mtpCoefs = MetaCyclicsM(m, maxOrder);
        foreach (var (_, coefs) in mtpCoefs)
        {
            foreach (var (_, n, r) in coefs)
            {
                ++total;
                var mt = FG.MetaCyclicPg(m, n, r);
                var (_, abType, m2, n2, r2) = RewriteMetaCyclic(m, n, r);
                var mt2 = FG.MetaCyclicPg(m2, n2, r2);
                var ab = FG.AbelianPerm(abType);
                var prod = ab.Order == 1 ? mt : ProductPermGroup(ab, mt2);
                DisplayGroup.HeadOrders(mt);
                var subs = mt.AllSubgroups();
                var decomp = subs.DecomposeProducts(subs.ProperNonTrivialNormalSubgroups());
                if (ab.Order != 1)
                {
                    DisplayGroup.HeadOrders(prod);
                    ++found;
                    if (!GroupCraft.AreIsomorphic(mt, prod))
                        throw new();

                    var (lhs, rhs, _) = decomp.OrderByDescending(e => e.lhs.Order)
                        .FirstOrDefault(e =>
                                e.isDirectProduct && (e.lhs.GroupType == GroupType.AbelianGroup ||
                                                      e.rhs.GroupType == GroupType.AbelianGroup),
                            (lhs: subs.TrivialClass, rhs: subs.TrivialClass, isDirectProduct: true));

                    if (lhs.GroupType == GroupType.NonAbelianGroup)
                        (lhs, rhs) = (rhs, lhs);
                    if (lhs.Order == 1)
                    {
                        subs.Naming();
                        decomp.OrderByDescending(e => e.lhs.Order).Where(e => e.isDirectProduct).Take(5).Println();
                        throw new();
                    }

                    var abName = Group.AbelianGroupType(lhs.Representative).ToAbString();
                    var sdp = subs.Restriction(rhs.Representative);
                    sdp.Naming();
                    Console.WriteLine($"{mt} = {abName} x {sdp.Parent}");
                    Console.WriteLine();
                }
                else
                {
                    var (lhs, rhs, _) = decomp.OrderByDescending(e => e.lhs.Order)
                        .FirstOrDefault(e =>
                                e.isDirectProduct && (e.lhs.GroupType == GroupType.AbelianGroup ||
                                                      e.rhs.GroupType == GroupType.AbelianGroup),
                            (lhs: subs.TrivialClass, rhs: subs.TrivialClass, isDirectProduct: true));
                    if (lhs.Order != 1)
                    {
                        if (lhs.GroupType == GroupType.NonAbelianGroup)
                            (lhs, rhs) = (rhs, lhs);
                        var abName = Group.AbelianGroupType(lhs.Representative).ToAbString();
                        var sdp = subs.Restriction(rhs.Representative);
                        sdp.Naming();
                        Console.WriteLine($"{mt} = {abName} x {sdp.Parent}");
                        throw new();
                    }
                    else
                    {
                        ++notFound;
                        Console.WriteLine($"M x M(m',n')r' not found for {mt}");
                        Console.WriteLine();
                    }
                }
            }
        }
    }

    Console.WriteLine($"Total:{total} Found:{found} NotFound:{notFound}");
    Console.WriteLine();
    if (total != found + notFound)
        throw new();
}

void DecomposeMetaCyclic(int maxM, int maxOrder)
{
    var (total, found) = (0, 0);
    foreach (var m in (maxM - 2).SeqLazy(3))
    {
        var mtpCoefs = MetaCyclicsM(m, maxOrder);
        foreach (var (_, coefs) in mtpCoefs)
        {
            foreach (var (_, n, r) in coefs)
            {
                ++total;
                var (_, abType, m2, n2, r2) = RewriteMetaCyclic(m, n, r);
                var nameMt = FG.MetaCyclicName(m, n, r);
                var nameMt2 = FG.MetaCyclicName(m2, n2, r2);
                var nameAb = abType.ToAbString();
                var ordAb = abType.Aggregate(1, (acc, e) => e * acc);
                if (ordAb != 1)
                {
                    ++found;
                    Console.WriteLine($"{nameMt,-15} = {nameMt2,-15} x {nameAb}");
                }
                else
                    Console.WriteLine(nameMt);
            }
        }
    }

    Console.WriteLine($"Total:{total} Found:{found}");
    Console.WriteLine();
}

void RunDecomposeMetacyclic()
{
    TestDecomposeMetaCyclic(maxM: 32, maxOrder: 128); // Total:302 Found:178 NotFound:124
    DecomposeMetaCyclic(maxM: 64, maxOrder: 1024); // Total:4121 Found:3375

    var ord = 5 * 9 * 4;
    FG.MetaCyclicCoefs(ord).Where(e => e.r > 1)
        .ToDictionary(e => (e.m, e.n, e.r, details: FG.MetaCyclicGens(e.m, e.n, e.r)),
            e => RewriteMetaCyclic(e.m, e.n, e.r).name)
        .Println(e => $"{e.Key.details.name,-20} = {e.Value,-30} gens:{e.Key.details.gens.ToXSet()}",
            $"All Metacyclic groups of order {ord}");
    // All Metacyclic groups of order 180
    // M(3x:60)2            = C15 x Dic3                     gens:{ [(2 3)(4 5 6)(7 8 9 10)(11 12 13 14 15)], [(1 2 3)] }
    // F(5x:36)2            = C9 x F(5x:4)2                  gens:{ [(2 4 5 3)(6 7 8 9 10 11 12 13 14)], [(1 2 3 4 5)] }
    // F(5x:36)4            = C9 x Dic5                      gens:{ [(2 5)(3 4)(6 7 8 9)(10 11 12 13 14 15 16 17 18)], [(1 2 3 4 5)] }
    // M(6x:30)5            = C30 x D6                       gens:{ [(4 5)(6 7 8)(9 10 11 12 13)], [(1 2)(3 4 5)] }
    // F(9x:20)8            = C5 x Dic9                      gens:{ [(2 9)(3 8)(4 7)(5 6)(10 11 12 13)(14 15 16 17 18)], [(1 2 3 4 5 6 7 8 9)] }
    // M(10x:18)9           = C18 x D10                      gens:{ [(4 7)(5 6)(8 9 10 11 12 13 14 15 16)], [(1 2)(3 4 5 6 7)] }
    // M(15x:12)2           = C3 x F(15x:4)2                 gens:{ [(2 3)(5 7 8 6)(9 10 11)], [(1 2 3)(4 5 6 7 8)] }
    // M(15x:12)4           = C3 x C3 x Dic5                 gens:{ [(5 8)(6 7)(9 10 11)(12 13 14 15)], [(1 2 3)(4 5 6 7 8)] }
    // M(15x:12)7           = C3 x C3 x F(5x:4)2             gens:{ [(5 7 8 6)(9 10 11)], [(1 2 3)(4 5 6 7 8)] }
    // M(15x:12)11          = C15 x Dic3                     gens:{ [(2 3)(9 10 11)(12 13 14 15)], [(1 2 3)(4 5 6 7 8)] }
    // M(15x:12)14          = C3 x Dic15                     gens:{ [(2 3)(5 8)(6 7)(9 10 11)(12 13 14 15)], [(1 2 3)(4 5 6 7 8)] }
    // M(18x:10)17          = C10 x D18                      gens:{ [(4 11)(5 10)(6 9)(7 8)(12 13 14 15 16)], [(1 2)(3 4 5 6 7 8 9 10 11)] }
    // M(30x:6)11           = C30 x D6                       gens:{ [(4 5)(11 12 13)], [(1 2)(3 4 5)(6 7 8 9 10)] }
    // M(30x:6)19           = C6 x C3 x D10                  gens:{ [(7 10)(8 9)(11 12 13)], [(1 2)(3 4 5)(6 7 8 9 10)] }
    // M(30x:6)29           = C6 x D30                       gens:{ [(4 5)(7 10)(8 9)(11 12 13)], [(1 2)(3 4 5)(6 7 8 9 10)] }
    // F(45x:4)8            = F(45x:4)8                      gens:{ [(2 3 5 4)(7 14)(8 13)(9 12)(10 11)], [(1 2 3 4 5)(6 7 8 9 10 11 12 13 14)] }
    // M(45x:4)19           = C9 x Dic5                      gens:{ [(2 5)(3 4)(15 16 17 18)], [(1 2 3 4 5)(6 7 8 9 10 11 12 13 14)] }
    // M(45x:4)26           = C5 x Dic9                      gens:{ [(7 14)(8 13)(9 12)(10 11)(15 16 17 18)], [(1 2 3 4 5)(6 7 8 9 10 11 12 13 14)] }
    // M(45x:4)28           = C9 x F(5x:4)3                  gens:{ [(2 3 5 4)], [(1 2 3 4 5)(6 7 8 9 10 11 12 13 14)] }
    // Dic45                = Dic45                          gens:{ [(2 5)(3 4)(7 14)(8 13)(9 12)(10 11)(15 16 17 18)], [(1 2 3 4 5)(6 7 8 9 10 11 12 13 14)] }
    // M(90x:2)19           = C18 x D10                      gens:{ [(4 7)(5 6)], [(1 2)(3 4 5 6 7)(8 9 10 11 12 13 14 15 16)] }
    // M(90x:2)71           = C10 x D18                      gens:{ [(9 16)(10 15)(11 14)(12 13)], [(1 2)(3 4 5 6 7)(8 9 10 11 12 13 14 15 16)] }
    // D180                 = C2 x D90                       gens:{ [(4 7)(5 6)(9 16)(10 15)(11 14)(12 13)], [(1 2)(3 4 5 6 7)(8 9 10 11 12 13 14 15 16)] }
}

void testMt()
{
    var (total, found) = (0, 0);
    for (int k = 1; k <= 8; k++)
    {
        ++total;
        var (m, n, r) = (15, 4 * k, 8);
        var mt = FG.MetaCyclicPg(m, n, r);
        DisplayGroup.HeadOrdersGenerators(mt);
        var autMt = Group.AutomorphismGroup(mt);
        DisplayGroup.HeadOrders(autMt);

        var (name2, abType, m1, n1, r1) = RewriteMetaCyclic(m, n, r);
        Console.WriteLine($"{mt} => {name2}");
        var (name1, sn1, gens1) = FG.MetaCyclicGens(m1, n1, r1);
        var mt1 = Group.Generate(name1, sn1, gens1);
        var autMt1 = RegPermAutGroup(mt1);
        var innAutMt1 = RegPermInnAutGroup(mt1);
        DisplayGroup.HeadOrdersGenerators(mt1);
        DisplayGroup.HeadOrders(autMt1);
        DisplayGroup.HeadNames(innAutMt1);
        DisplayGroup.Generators(innAutMt1);

        var (autMt2, map) = autMt1.ToPermGroup();
        var innAutMt2 = Group.Generate(innAutMt1.Name, autMt2.BaseGroup,
            innAutMt1.GetGenerators().Select(e => map[e]).ToArray());
        foreach (var (K, conj) in GroupCraft.AllFactors(autMt2, innAutMt2)
                     .DistinctBy(e => e.Key.ElementsOrdersList().Glue(",")))
        {
            DisplayGroup.HeadOrdersNames(K.ToPermGroup().Item1);
            ++found;
            break;
        }

        Console.WriteLine();
    }

    Console.WriteLine(new { total, found });
}

ConcreteGroup<Perm> HolCp(int p)
{
    if (IntExt.PrimesDec(p).Count != 1)
        throw new();

    var a = FG.Cycles(p);
    var upGens = FG.UnPermGens(p);
    return Group.Generate($"Hol[C{p}]", a.Sn, [a, ..upGens]);
}

ConcreteGroup<Perm> HolCn(int n)
{
    var a = FG.Cycles(n);
    return Group.Generate($"Hol[C{n}]", a.Sn, [a, ..FG.UnPermGens(n)]);
    // var pks = IntExt.PrimesDec(n).Select(e => e.Key.Pow(e.Value)).Order().ToArray();
    // return ProductPermGroup(pks.Select(e => HolCp(e)).ToArray());
}

(string name, ConcreteGroup<Perm> mt, ConcreteGroup<Perm> autMt) AutFrob(int m, int n, int r)
{
    var mt = FG.MetaCyclicPg(m, n, r);
    var (name, abType, m1, n1, r1) = RewriteMetaCyclic(m, n, r);
    var mt1 = FG.MetaCyclicPg(m1, n1, r1);
    var ab = FG.AbelianPerm(abType);
    if (IntExt.Gcd(m1, n1) != 1 || IntExt.Gcd(ab.Order, n1) != 1)
        throw new();

    var autAb = RegPermAutGroup(ab);
    var z1Type = IntExt.Gcd(m, n) == 1 ? Group.AbelianGroupType(RegPermZentrumAutGroup(mt1)) : [1];
    var z1pg = FG.AbelianPerm(z1Type);
    var head = ProductPermGroup(autAb, z1pg);
    if (head.GroupType == GroupType.AbelianGroup)
        head = FG.AbelianPerm(Group.AbelianGroupType(head));

    var tail = HolCn(m1);
    var autMt = head.Order != 1 ? ProductPermGroup([head, tail]) : tail;
    return (name, mt, autMt);
}

IEnumerable<(string name, ConcreteGroup<Perm> mt, ConcreteGroup<Perm> autMt)> RunAutFrobM(int m, int maxOrd)
{
    foreach (var (_, coefs) in MetaCyclicsM(m, maxOrd))
    foreach (var (_, n, r) in coefs.Where(e => IntExt.Gcd(m, e.n) == 1))
        yield return AutFrob(m, n, r);
}

bool AutFrobVerification((string name, ConcreteGroup<Perm> mt, ConcreteGroup<Perm> autMt) fAut)
{
    var (name, mt, autMt) = fAut;
    Console.WriteLine($"##################################### {mt} = {name}");
    DisplayGroup.HeadOrders(mt);
    DisplayGroup.HeadOrders(autMt);

    var autMt2 = RegPermAutGroup(mt);
    DisplayGroup.HeadOrders(autMt2);
    Console.WriteLine();
    var autMtGens = autMt.GetGenerators().Count();
    var autMt2Gens = autMt2.GetGenerators().Count();
    if (autMtGens <= 3 || autMt2Gens <= 3)
    {
        if (autMt2Gens <= 3 && GroupCraft.AreIsomorphic(autMt2, autMt))
            return true;

        if (autMtGens <= 3 && GroupCraft.AreIsomorphic(autMt, autMt2))
            return true;
    }
    else if (autMt.ElementsOrdersList().SequenceEqual(autMt2.ElementsOrdersList()))
        return true;

    return false;
}

void RunAutFrobWithVerification(int maxM, int maxOrd)
{
    GlobalStopWatch.Restart();

    var (total, found) = (0, 0);
    foreach (var m in (maxM - 2).SeqLazy(3))
    {
        foreach (var fAut in RunAutFrobM(m, maxOrd))
        {
            ++total;
            var check = AutFrobVerification(fAut);
            if (check)
                ++found;
        }
    }

    GlobalStopWatch.Show($"Total:{total} AutFound:{found}");
}

void RunAutFrob(int maxM, int maxOrd)
{
    GlobalStopWatch.Restart();

    var total = 0;
    foreach (var m in (maxM - 2).SeqLazy(3))
    {
        foreach (var (name, mt, autMt) in RunAutFrobM(m, maxOrd))
        {
            ++total;
            Console.WriteLine($"##################################### {mt} = {name}");
            DisplayGroup.HeadOrders(mt);
            DisplayGroup.HeadOrders(autMt);
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show($"Total:{total}");
}

void DisplayAutFrob(int m, int n, int r)
{
    var (name, mt, autMt) = AutFrob(m, n, r);
    Console.WriteLine($"##################################### {mt} = {name}");
    DisplayGroup.HeadOrders(mt);
    DisplayGroup.HeadOrders(autMt);
    Console.WriteLine();
}

(string name, Perm a, Perm b, Dictionary<Perm, (int i, int j)>) DecomposeMetacyclic(int m, int n, int r)
{
    var (name, sn, gens) = FG.MetaCyclicGens(m, n, r);
    var (a, b) = (gens[0], gens[1]);
    var table = m.Range().Grid2D(n.Range()).ToDictionary(e => (a ^ e.t1) * (b ^ e.t2), e => e);
    return (name, a, b, table);
}

{
    // M(m,n)r = Ab x M(m',n')r'
    // with m and n coprimes
    // m' = p1^k1 x ... x pl^kl
    // Aut[M(m,n)r] = Aut[Ab] x Z(Aut[M(m',n')r']) x Prod[Hol[Cpi^ki]]
    // M(21x:8)13 = C3 x F(7x:8)6 and Aut[M(21x:8)13] = C2 x C2 x C2 x Hol[C7]
    // M(15x:16)13 = C3 x F(5x:16)3 and Aut[M(15x:16)13] = C2 x C4 x Hol[C5]
    // RunAutFrobWithVerification(maxM: 32, maxOrd: 256); // Total:235 AutFound:235 Time:1m22s

    // RunAutFrob(maxM: 64, maxOrd: 128); // # Total:151 Time:14.415s
}

void AutMetacyclic(int maxM, int maxOrd)
{
    GlobalStopWatch.Restart();
    var (total, found, subTotal1, subTotal2) = (0, 0, 0, 0);
    foreach (var m in (maxM - 2).SeqLazy(3))
    {
        foreach (var (r, coefs) in MetaCyclicsM(m, maxOrd))
        {
            foreach (var (_, n, _) in coefs)
            {
                var (nameAbMt1, abType, m1, n1, r1) = RewriteMetaCyclic(m, n, r);
                if (IntExt.Gcd(m1, n1) != 1)
                    continue;

                ++total;

                var mt = FG.MetaCyclicPg(m, n, r);
                var autMt = RegPermAutGroup(mt);
                var holCm1 = HolCn(m1);
                
                var iso = GroupCraft.AllIsomorphicSubGroups(autMt, holCm1).First(iso => GroupCraft.IsNormal(autMt, iso));
                var quo = autMt.Over(iso);
                var Quo = GroupCraft.DummyFactor(autMt, iso);
                if (quo.Order == Quo.Order && iso.Intersect(Quo).Count() == 1)
                {
                    ++found;

                    var ab = FG.AbelianPerm(abType);
                    var autAb = RegPermAutGroup(ab);
                    if (Quo.Order % autAb.Order == 0)
                    {
                        Console.WriteLine($"##################################### {mt.Name} = {nameAbMt1}");
                        DisplayGroup.HeadOrders(mt);
                        DisplayGroup.HeadOrders(autMt);
                        DisplayGroup.HeadOrders(iso);
                        if (Quo.Order == 1)
                            continue;
                        
                        DisplayGroup.HeadOrders(Quo);
                        GroupCraft.DisplayIsNormal(autMt, iso);
                        GroupCraft.DisplayIsNormal(autMt, Quo);
                        if (Quo.GroupType == GroupType.AbelianGroup || autAb.Order == 1)
                            continue;
                        
                        var (iso2, Quo2) = GroupCraft.AllIsomorphicSubGroups(Quo, autAb)
                            .Select(iso2 => (iso2, quo2: GroupCraft.DummyNormal(Quo, iso2)))
                            .First(e => Quo.Order == e.iso2.Order * e.quo2.Order && e.iso2.Intersect(e.quo2).Count() == 1);
                        DisplayGroup.HeadOrders(Quo2);
                        DisplayGroup.HeadOrders(iso2);
                        GroupCraft.DisplayIsNormal(Quo, Quo2);
                        GroupCraft.DisplayIsNormal(Quo, iso2);
                        if (iso2.GroupType == GroupType.NonAbelianGroup || Quo2.GroupType == GroupType.NonAbelianGroup)
                        {
                            ++subTotal1;
                            Console.WriteLine("warning");
                        }
                        else
                            ++subTotal2;
                        
                        Console.WriteLine();
                    }
                }
            }
        }
    }

    GlobalStopWatch.Show($"Total:{total} SubTotal:{subTotal1}/{subTotal2} Found:{found}");
}

{
    // AutMetacyclic(32, 128); // Total:192
    AutMetacyclic(64, 128);
    // AutMetacyclic(32, 256);
}