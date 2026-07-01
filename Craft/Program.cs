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
    var bgAut = new AutomorphismGroup<T>(g);
    var allAut = GroupCraft.AllMorphismsWithPruning(g, g, Group.MorphismType.Isomorphism)
        .Select(aut => new Automorphism<T>(bgAut, aut.HomMap))
        .OrderByDescending(aut => Group.GenerateElements(bgAut, aut).Count)
        .ToArray();
    var autG = Group.Generate($"Aut[{g.Name}]", bgAut, allAut);
    return autG;
}

ConcreteGroup<Perm> RegPermAutGroup<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
{
    var bgAut = new AutomorphismGroup<T>(G);
    var sn = new Sn(G.Count());
    var mapEltIdx = G.OrderBy(e => G.ElementsOrders[e]).Index().ToDictionary(e => e.Item, e => e.Index);
    var allAut = GroupCraft.AllMorphismsWithPruning(G, G, Group.MorphismType.Isomorphism)
        .Select(aut => new Automorphism<T>(bgAut, aut.HomMap))
        .Select(e => e.AutMap.ToDictionary(f => mapEltIdx[f.Key], f => mapEltIdx[f.Value]))
        .Select(e => sn.CreateElementTable(e.OrderBy(f => f.Key).Select(f => f.Value).ToArray()))
        .ToArray();
    var autG = Group.Generate($"Aut[{G.Name}]", sn, allAut);
    return autG;
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
    var autG = Group.Generate($"InnAut[{G.Name}]", sn, gensInnAut);
    return autG;
}

(ConcreteGroup<Perm> G, ConcreteGroup<Perm> AutG) RegPermGroupAndAutGroup<T>(ConcreteGroup<T> G)
    where T : struct, IElt<T>
{
    var sn = new Sn(G.Count());
    var mapEltIdx = G.Index().OrderBy(e => G.ElementsOrders[e.Item])
        .ToDictionary(e => e.Item, e => e.Index);
    var mapIdxElt = mapEltIdx.ToDictionary(e => e.Value, e => e.Key);
    GroupAction<T, Perm> act = (e, p) =>
        sn.CreateElementTable(p.Table.Select(i => mapEltIdx[G.Op(e, mapIdxElt[i])]).ToArray());
    var mapReg = G.ToDictionary(e => e, e => act(e, sn.Neutral()));
    var gensReg = G.GetGenerators().Select(e => mapReg[e]).ToArray();
    var Greg = Group.Generate(G.Name, sn, gensReg);

    var bgAut = new AutomorphismGroup<T>(G);
    var allAut = GroupCraft.AllMorphismsWithPruning(G, G, Group.MorphismType.Isomorphism)
        .Select(aut => new Automorphism<T>(bgAut, aut.HomMap))
        .Select(e => e.AutMap.ToDictionary(f => mapEltIdx[f.Key], f => mapEltIdx[f.Value]))
        .Select(e => sn.CreateElementTable(e.OrderBy(f => f.Key).Select(f => f.Value).ToArray()))
        .OrderByDescending(e => e.Order)
        .ToArray();
    var autG = Group.Generate($"Aut[{G.Name}]", sn, allAut);
    return (Greg, autG);
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

void AutMMnames()
{
    for (int k = 4; k <= 8; k++)
    {
        var mm = FG.ModularMaxPg(k);
        var autMM = RegPermGroupAndAutGroup(mm).AutG;
        DisplayGroup.HeadOrders(autMM);
        DisplaySDPAb(autMM.AllSubgroups());
    }
}

void AutQDnames()
{
    for (int k = 4; k <= 7; k++)
    {
        var mm = FG.SemiDihedralPg(k);
        var autQD = RegPermGroupAndAutGroup(mm).AutG;
        DisplayGroup.HeadOrders(autQD);
        DisplaySDPAb(autQD.AllSubgroups());
    }
}

void AutMMSDP()
{
    for (int k = 4; k <= 8; k++)
    {
        var n = 2.Pow(k);
        var mm = FG.ModularMaxPg(k);
        var sn = mm.Neutral().Sn;
        DisplayGroup.HeadOrders(mm);
        var autMM = RegPermAutGroup(mm);
        DisplayGroup.HeadOrders(autMM);
        var l = n / 8;
        var (clc2c2, autclc2c2) = RegPermGroupAndAutGroup(FG.AbelianPerm(l, 2, 2));
        var c2c2 = FG.AbelianPerm(2, 2);
        var counter = 0;
        foreach (var theta in GroupCraft.AllMorphismsWithPruning(c2c2, autclc2c2))
        {
            ++counter;
            var gensAutMM = clc2c2.GetGenerators().Concat(c2c2.GetGenerators().Select(e => theta[e])).ToArray();
            var autMMpg = Group.Generate(autMM.Name, sn, gensAutMM);
            if (autMMpg.ElementsOrdersList().SequenceEqual(autMM.ElementsOrdersList()))
            {
                DisplayGroup.HeadOrdersGenerators(autMMpg);
                break;
            }
        }

        Console.WriteLine($"try:{counter}");
        Console.WriteLine();
    }
}

void AutQDSDP()
{
    for (int k = 4; k <= 7; k++)
    {
        var n = 2.Pow(k);
        var qd = FG.SemiDihedralPg(k);
        DisplayGroup.HeadOrders(qd);
        var autQD = RegPermAutGroup(qd);
        DisplayGroup.HeadOrders(autQD);
        var l1 = n / 4;
        var l2 = n / 8;
        var (cl, autcl) = RegPermGroupAndAutGroup(FG.AbelianPerm(l1));
        var cl2c2 = FG.AbelianPerm(l2, 2);
        var sn1 = cl.Neutral().Sn;
        var sn2 = cl2c2.Neutral().Sn;
        var counter = 0;
        foreach (var theta in GroupCraft.AllMorphismsWithPruning(cl2c2, autcl))
        {
            ++counter;
            var gensLhs = cl.GetGenerators().Select(e => FG.PaddingRight(e, sn2.N)).ToArray();
            var gensRhs = cl2c2.GetGenerators().Select(e => FG.ConcatPerm(theta[e], e)).ToArray();
            var gensAutQD = gensLhs.Concat(gensRhs).ToArray();
            var autQDpg = Group.Generate(autQD.Name, gensAutQD[0].Sn, gensAutQD);
            if (autQDpg.ElementsOrdersList().SequenceEqual(autQD.ElementsOrdersList()))
            {
                DisplayGroup.HeadOrdersGenerators(autQDpg);
                break;
            }
        }

        Console.WriteLine($"try:{counter}");
        Console.WriteLine();
    }
}

void RunAutMM()
{
    for (int k = 4; k <= 10; k++)
    {
        var n = 2.Pow(k);
        var mm = FG.ModularMaxPg(k);
        DisplayGroup.HeadOrders(mm);
        var l = n / 8;
        var sm = mm.Neutral().Sn;
        var m = sm.N;

        var a = sm.OpSeq(4.SeqLazy().Select(i => sm.CycleP1((m / 4).SeqLazy(i, 4).ToArray())));
        var b = sm.OpSeq((m / 2).SeqLazy(0, 2).Select(i => sm.CycleP1([i, i + 1])).ToArray());
        var c = sm.OpSeq(2.SeqLazy().SelectMany(i => (m / 4).SeqLazy(i, 4).Select(j => sm.CycleP1([j, j + 2])))
            .ToArray());
        var d = sm.OpSeq(a.DisjoinCycles.Take(2).Select(e => e ^ (m / 8)));

        var sdp = Group.Generate($"(C{l} x C2 x C2) x: C2", sm, a, b, c, d);
        DisplayGroup.HeadOrdersGenerators(sdp);
        var autMM = RegPermAutGroup(mm);
        DisplayGroup.HeadOrders(autMM);
    }
}

void RunAutQD()
{
    for (int k = 4; k <= 8; k++)
    {
        var n = 2.Pow(k);
        var qd = FG.SemiDihedralPg(k);
        DisplayGroup.HeadOrders(qd);
        var l1 = n / 4;
        var l2 = n / 8;
        var (a0, b0, c0) = (FG.Cycles(l1), FG.Cycles(l2), FG.Cycles(2));
        var sn1 = a0.Sn;
        var sn2 = b0.Sn;
        var sn3 = c0.Sn;

        var a = FG.PaddingRight(FG.Cycles(l1), sn2.N + sn3.N);
        var b = FG.Padding(sn1.N, b0, sn3.N);
        var c1 = sn1.OpSeq((l1 / 2 - 1).SeqLazy(1).Select(i => sn1.CycleP1([i, l1 - i])));
        var c = FG.ConcatPerm(FG.PaddingRight(c1, sn2.N), c0);
        var sdp = Group.Generate($"C{l1} x: (C{l2} x C2)", a.Sn, a, b, c);
        DisplayGroup.HeadOrders(sdp);
        var autQD = RegPermAutGroup(qd);
        DisplayGroup.HeadOrders(autQD);
    }
}

void TestMM()
{
    AutMMnames();
    AutMMSDP();
    RunAutMM();
}

void TestQD()
{
    AutQDnames();
    AutQDSDP();
    RunAutQD();
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

// {
//     MetaCyclicsM(15, 512).Values.SelectMany(e => e)
//         .Select(e => $"{FG.MetaCyclicName(e.m, e.n, e.r),-20} = {RewriteMetaCyclic(e.m, e.n, e.r).name}")
//         .Println();
// }

ConcreteGroup<Perm> HolCp(int p)
{
    if (IntExt.PrimesDec(p).Count != 1)
        throw new();

    var q = IntExt.Phi(p);
    var rq = IntExt.Coprimes(p).OrderDescending().First(rq => rUmToPerm(p, rq).Order == q);
    var mtpq = FG.MetaCyclicPg(p, q, rq);
    mtpq.Name = $"Hol[C{p}]";
    return mtpq;
}

void runInnOutAut(int maxM, int maxOrd)
{
    GlobalStopWatch.Restart();

    var (total, found) = (0, 0);
    foreach (var m in (maxM - 2).SeqLazy(3))
    {
        foreach (var (_, coefs) in MetaCyclicsM(m, maxOrd))
        {
            foreach (var (_, n, r) in coefs.Where(e => IntExt.Gcd(e.m, e.n) == 1))
            {
                ++total;
                Console.WriteLine($"########################################## Start m:{m} n:{n} r:{r}");
                var mt = FG.MetaCyclicPg(m, n, r);
                var autMt = RegPermAutGroup(mt);
                var (name, abType, m1, n1, r1) = RewriteMetaCyclic(m, n, r);
                var mt1 = FG.MetaCyclicPg(m1, n1, r1);
                var autMt1 = RegPermAutGroup(mt1);

                var autAb = RegPermAutGroup(FG.AbelianPerm(abType));
                if (abType.Length == 1)
                    autAb.Name = Group.AbelianGroupType(autAb).ToAbString();
                
                var decpk = IntExt.PrimesDec(m1).Select(e => e.Key.Pow(e.Value)).Order().ToArray();
                var zType = Group.AbelianGroupType(Group.Zentrum(autMt1));
                var zpg = FG.AbelianPerm(zType);
                var autMt2 = ProductPermGroup(decpk.Select(e => HolCp(e)).Prepend(zpg).Prepend(autAb)
                    .Where(e => e.Order != 1).ToArray());

                Console.WriteLine($"{mt} = {name}");
                DisplayGroup.HeadOrders(mt);
                DisplayGroup.HeadOrders(autMt);
                DisplayGroup.HeadOrders(autMt2);
                if (mt.Order <= maxOrd / 2)
                {
                    if (GroupCraft.AreIsomorphic(autMt, autMt2))
                        ++found;
                }
                else if (autMt.ElementsOrdersList().SequenceEqual(autMt2.ElementsOrdersList()))
                    ++found;

                Console.WriteLine();
            }
        }
    }

    GlobalStopWatch.Show($"Total:{total} AutFound:{found}");
}

{
    // M(m,n)r = Ab x M(m',n')r'
    // with m and n coprimes
    // m' = p1^k1 x ... x pl^kl
    // Aut[M(m,n)r] = Aut[Ab] x Z(Aut[M(m',n')r']) x Prod[Hol[Cpi^ki]]
    // M(21x:8)13 = C3 x F(7x:8)6 and Aut[M(21x:8)13] = C2 x C2 x C2 x Hol[C7]
    // M(15x:16)13 = C3 x F(5x:16)3 and Aut[M(15x:16)13] = C2 x C4 x Hol[C5]
    runInnOutAut(maxM: 32, maxOrd: 256); // Total:235 AutFound:235 Time:1m22s
}