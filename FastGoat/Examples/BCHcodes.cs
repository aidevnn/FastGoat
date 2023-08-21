using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

/// <summary>
/// Provides methods for working with BCH (Bose-Chaudhuri-Hocquenghem) codes,
/// including encoding and decoding using the Berlekamp-Massey algorithm.
/// Algebre Tome1, page 387-397
/// </summary>
public static class BCHcodes
{
    /// <summary>
    /// Calculates the cyclotomic classes for a given prime modulus and order.
    /// </summary>
    /// <param name="n">The order of the abelian group.</param>
    /// <param name="q">Cl(i)~Cl(j) iff j = i*q^k .</param>
    /// <returns>An array of hash sets representing the cyclotomic classes.</returns>
    static HashSet<int>[] CyclotomicClasses(int q, int n)
    {
        var cn = new Cn(n);
        var set = new HashSet<HashSet<int>>(new SetEquality<int>());
        var q0 = new ZnInt(n, q);
        foreach (var i0 in cn)
        {
            var Si = new HashSet<int>() { i0.K };
            var s = 1;
            while (true)
            {
                var qsi = q0.Pow(s) * i0;
                if (qsi.Equals(i0))
                    break;

                Si.Add(qsi.K);
                ++s;
            }

            set.Add(Si.ToHashSet());
        }

        return set.ToArray();
    }

    /// <summary>
    /// Generates a tuple containing the Galois field Fq of order n=q-1
    /// and the generator polynomial KPoly&lt;ZnInt&gt; for the BCH codes
    /// </summary>
    /// <param name="m">The parameter 'm' representing the exponent of the Galois field, Fq with q=2^m.</param>
    /// <param name="delta">The parameter 'delta' representing the Hamming distance.</param>
    /// <returns>A tuple containing the Fq field and the KPoly&lt;ZnInt&gt; polynomial generator.</returns>
    static (Fq fq, KPoly<ZnInt> g) GeneratorBCHmδ(int m, int delta)
    {
        var fq = new Fq(2.Pow(m), 'a');
        var a = fq.X;
        var X = FG.KPoly('X', fq.X);
        // var nth = Group.MulGroup(fq.Name, a);
        // DisplayGroup.HeadElements(nth);
        var n = fq.Q - 1;
        var q = fq.P;

        var set = CyclotomicClasses(q, n);
        var set1 = set.Select(si => (si, si.Select(k => X - a.Pow(k)).Aggregate((s0, s1) => s0 * s1))).OrderBy(si => si.Item2)
            .ToArray();
        set1.Println(si => $"[{si.si.Glue("; ")}] => {si.Item2}", $"{fq.FullName} nb cyclo classes = {set.Length}");

        var rg = (delta - 1).Range(1).Select(e0 => set1.First(e1 => e1.si.Contains(e0)).Item2).ToArray();
        var codeBCH = rg.Aggregate((e0, e1) => e0 * e1 / Ring.Gcd(e0, e1));
        rg.Println($"BCH(q:{q}, n:{n}, δ:{delta}) Code = {codeBCH}");
        Console.WriteLine();

        if (codeBCH.Coefs.Any(c => c.Poly.Degree != 0))
            throw new();

        return (fq, new('x', ZnInt.ZnZero(q), codeBCH.Coefs.Select(c => c[0]).ToArray()));
    }

    /// <summary>
    /// Generates an array of random binaries (ZnInt modulus 2) values representing a word of given length.
    /// </summary>
    /// <param name="k">The length of the word.</param>
    /// <returns>An array of binaries values representing the random word.</returns>
    static ZnInt[] RandWord(int k) => k.Range().Select(i => new ZnInt(2, IntExt.Rng.Next(2))).TrimSeq().ToArray();

    /// <summary>
    /// Converts an array of binaries (ZnInt modulus 2) values representing a word into a KPoly&lt;ZnInt&gt; polynomial.
    /// </summary>
    /// <param name="word">The array of binaries values representing the word.</param>
    /// <returns>A KPoly&lt;ZnInt&gt; polynomial representing the word.</returns>
    static KPoly<ZnInt> WordPoly(ZnInt[] word) => new('x', ZnInt.ZnZero(2), word);

    /// <summary>
    /// Encodes a given polynomial word using the BCH (Bose-Chaudhuri-Hocquenghem) encoding technique.
    /// </summary>
    /// <param name="word">The polynomial representing the input word with binary values.</param>
    /// <param name="g">The generator polynomial used for BCH codes.</param>
    /// <returns>A KPoly&lt;ZnInt&gt; polynomial representing the encoded word.</returns>
    static KPoly<ZnInt> EncodeBCH(KPoly<ZnInt> word, KPoly<ZnInt> g) => word * g;

    /// <summary>
    /// Converts a polynomial with binary values to a binary string representation, padded to a specified word length.
    /// </summary>
    /// <param name="P">The polynomial with binary values to be converted.</param>
    /// <param name="n">The desired length of the word.</param>
    /// <returns>A binary string representation of the polynomial, padded with '0' characters on the right
    /// to match the specified word length.</returns>
    static string Poly2Bin(KPoly<ZnInt> P, int n) => P.Coefs.Glue().PadRight(n, '0');

    /// <summary>
    /// Calculates the syndromes of a given polynomial 'm' over a Galois field Fq, for a specified Hamming distance 'delta'.
    /// </summary>
    /// <param name="m">The polynomial of binary values for which syndromes are to be calculated.</param>
    /// <param name="fq">The Galois field Fq with q = 2^m and order n = q - 1.</param>
    /// <param name="delta">The Hamming distance parameter.</param>
    /// <returns>A KPoly&lt;EPoly&lt;ZnInt&gt;&gt; polynomial representing the calculated syndromes.</returns>
    static KPoly<EPoly<ZnInt>> Syndrom(KPoly<ZnInt> m, Fq fq, int delta)
    {
        var b = fq.X;
        var X = FG.KPoly('X', b);
        return (delta - 1).Range(1).Select(j => m.Substitute(b.Pow(j)) * X.Pow(j)).Aggregate((e0, e1) => e0 + e1);
    }

    /// <summary>
    /// Introduces random noise into a given polynomial 'm' to simulate errors,
    /// with a specified Hamming distance 'delta'.
    /// </summary>
    /// <param name="m">The polynomial of binary values to which noise will be added.</param>
    /// <param name="delta">The Hamming distance parameter.</param>
    /// <param name="maxErrors">Indicates whether to generate the maximum number of errors
    /// or a random number of errors (default: false).</param>
    /// <returns>A KPoly&lt;ZnInt&gt; polynomial representing the original polynomial with added random noise.</returns>
    static KPoly<ZnInt> Noise(KPoly<ZnInt> m, int delta, bool maxErrors = false)
    {
        var d = m.Degree + 1;
        var t = (delta - 1) / 2;
        var t0 = maxErrors ? t - 1 : IntExt.Rng.Next(t);
        var m0 = new KPoly<ZnInt>(m.x, m.KZero, m.Coefs.ToArray());
        var errs = new HashSet<int>();
        while (errs.Count < t0)
        {
            var j = IntExt.Rng.Next(d);
            if (!errs.Add(j))
                continue;

            m0.Coefs[j] += 1;
        }

        return m0;
    }

    /// <summary>
    /// Applies the Berlekamp-Massey algorithm to decode a received noisy message
    /// by computing the localizator polynomial 'u' and the evaluator polynomial 'v'.
    /// </summary>
    /// <param name="S">The syndrom polynomial computed from the received noisy message.</param>
    /// <param name="delta">The Hamming distance parameter.</param>
    /// <returns>A tuple containing the localizator polynomial 'u' and the evaluator polynomial 'v'.</returns>
    static (KPoly<EPoly<ZnInt>> u, KPoly<EPoly<ZnInt>> v) BerlekampMassey(KPoly<EPoly<ZnInt>> S, int delta)
    {
        var t = (delta - 1) / 2;
        var (r0, u0, v0) = (S.One, S.Zero, S.X.Pow(2 * t));
        var (r1, u1, v1) = (S.Zero, S.One, S);

        // vi−1 = vi*qi + vi+1 puis les soustractions ri+1 = ri−1 − ri*qi et ui+1 = ui−1 − ui*qi
        // jusqu’à obtenir deg vi < t et deg vi−1 >= t
        while (!(v1.Degree < t && v0.Degree >= t))
        {
            var (quo, rem) = v0.Div(v1);
            (v0, v1) = (v1, rem);
            (r0, r1) = (r1, r0 - r1 * quo);
            (u0, u1) = (u1, u0 - u1 * quo);
        }

        var c = u1.LT;
        return (u1 / c, v1 / c);
    }

    /// <summary>
    /// Computes the error polynomial between the sent and received messages
    /// using the localizator polynomial 'u' and the Galois field 'Fq'.
    /// </summary>
    /// <param name="u">The localizator polynomial.</param>
    /// <param name="q">The Galois field Fq with order n=q-1.</param>
    /// <returns>A KPoly&lt;ZnInt&gt; polynomial representing the error between the sent and received messages.</returns>
    static KPoly<ZnInt> Error(KPoly<EPoly<ZnInt>> u, int q)
    {
        var b = u.KZero.X;
        var e = b.F.Zero.SubstituteChar('x');
        for (int i = 0; i < q - 1; i++)
        {
            var bi = b.Pow(i).Inv();
            if (u.Substitute(bi).IsZero())
                e += e.X.Pow(i);
        }

        return e;
    }

    /// <summary>
    /// Decodes a received polynomial with random noise using the Berlekamp-Massey algorithm and BCH codes.
    /// </summary>
    /// <param name="received">The received polynomial with random noise.</param>
    /// <param name="g">The generator polynomial of the BCH codes.</param>
    /// <param name="fq">The Galois Field Fq with q=2^m.</param>
    /// <param name="delta">The Hamming distance parameter.</param>
    /// <returns>A tuple containing the decoded word polynomial and the number of error bits.</returns>
    static (KPoly<ZnInt> word, int nbErr) DecodeBCH(KPoly<ZnInt> received, KPoly<ZnInt> g, Fq fq, int delta)
    {
        var S = Syndrom(received, fq, delta);
        var u = BerlekampMassey(S, delta).u;
        var err = Error(u, fq.Q);
        var decoded = received + err;
        var word = decoded / g;
        return (word, err.Coefs.Sum(e => e.K));
    }

    /// <summary>
    /// Runs multiple tests for BCH codes,
    /// including generating the generator polynomial, encoding words, adding noise, and decoding.
    /// </summary>
    /// <param name="m">The parameter 'm' for the Galois Field Fq with q=2^m.</param>
    /// <param name="delta">The Hamming distance parameter.</param>
    /// <param name="nbTest">The number of words to generate and test.</param>
    /// <param name="maxErrors">Indicates whether to generate the maximum number of errors
    /// or a random number of errors (default: false).</param>
    static void CheckBCHmδ(int m, int delta, int nbTest, bool maxErrors = false)
    {
        var (fq, g) = GeneratorBCHmδ(m, delta);
        var n0 = fq.Q - g.Degree - 1;
        Console.WriteLine($"BCH word {n0 + 1}bits, encoder {fq.Q}bits, max errors {(delta - 1) / 2 - 1}bits");
        Console.WriteLine($"code {g.Degree + 1}bits : {g}");
        Console.WriteLine();

        for (int k = 1; k <= nbTest; k++)
        {
            var word = WordPoly(RandWord(n0));
            var encoded = EncodeBCH(word, g);
            var received = Noise(encoded, delta, maxErrors);
            var (dword, err) = DecodeBCH(received, g, fq, delta);
            var verif = word.Equals(dword);
            Console.WriteLine($"Word{k,-3}    : {Poly2Bin(word, n0)}");
            Console.WriteLine($"  Send     : {Poly2Bin(encoded, fq.Q)}");
            Console.WriteLine($"  Received : {Poly2Bin(received, fq.Q)}");
            Console.WriteLine($"  Decoded  : {Poly2Bin(dword, n0)}");
            Console.WriteLine($"  NbErrors : {err}");
            Console.WriteLine($"  Verif    : {verif}");
            Console.WriteLine();
            if (!verif)
                throw new();
        }

        Console.WriteLine($"BCH word {n0 + 1}bits, encoder {fq.Q}bits, max errors {(delta - 1) / 2 - 1}bits");
        Console.WriteLine($"code {g.Degree + 1}bits : {g}");
        Console.WriteLine();
    }

    /// <summary>
    /// Runs multiple detailed tests for BCH codes
    /// including generating the generator polynomial, encoding words, adding noise, and decoding.
    /// </summary>
    /// <param name="m">The parameter 'm' for the Galois Field Fq with q=2^m.</param>
    /// <param name="delta">The Hamming distance parameter.</param>
    /// <param name="nbTest">The number of words to generate and test.</param>
    /// <param name="maxErrors">Indicates whether to generate the maximum number of errors
    /// or a random number of errors (default: false).</param>
    static void CheckDetailedBCHmδ(int m, int delta, int nbTest, bool maxErrors = false)
    {
        var (fq, g) = GeneratorBCHmδ(m, delta);
        var n0 = fq.Q - g.Degree - 1;
        Console.WriteLine($"BCH word {n0 + 1}bits, encoder {fq.Q}bits, max errors {(delta - 1) / 2 - 1}bits");
        Console.WriteLine($"code {g.Degree + 1}bits : {g}");
        Console.WriteLine();
        var t = (delta - 1) / 2;
        for (int k = 1; k <= nbTest; k++)
        {
            var send = EncodeBCH(WordPoly(RandWord(n0)), g);
            var received = Noise(send, delta, maxErrors);
            var err0 = received - send;
            var S = Syndrom(received, fq, delta);
            Console.WriteLine($"Syndrom : {S}");

            var b = S.KOne.X;
            var x = S.X;
            var I = err0.Coefs.Select((ei, i) => (ei, i)).Where(ei => !ei.ei.IsZero()).Select(ei => ei.i).ToArray();
            var u0 = I.Select(i => 1 - b.Pow(i) * x).Aggregate(x.One, (acc, xi) => acc * xi);
            var v0 = I.Select(i =>
                    b.Pow(i) * x * I.Except(new[] { i }).Select(j => 1 - b.Pow(j) * x).Aggregate(x.One, (acc, xj) => acc * xj))
                .Aggregate(x.Zero, (acc, xi) => acc + xi);

            var x2t = S.X.Pow(2 * t);
            var c = u0.LT;
            (u0, v0) = (u0 / c, v0 / c);
            Console.WriteLine($"  Direct          : u = {u0} and v = {v0}");
            var v2 = (u0 * S).Div(x2t).rem;
            if (!v0.Equals(v2))
                throw new();

            var (u, v) = BerlekampMassey(S, delta);
            Console.WriteLine($"  BerlekampMassey : u = {u} and v = {v}");

            var v1 = (u * S).Div(x2t).rem;
            Console.WriteLine($"  Verif Direct S  : {v0.Equals(v2)}");
            Console.WriteLine($"  Verif B-M    S  : {v.Equals(v1)}");
            if (!v.Equals(v1))
                throw new();
            Console.WriteLine($"  diff u : {u - u0}");
            Console.WriteLine($"  diff v : {v - v0}");

            var err1 = Error(u0, fq.Q);
            var err2 = Error(u, fq.Q);
            Console.WriteLine($"Send Word{k,-3}      : {Poly2Bin(send, fq.Q)}");
            Console.WriteLine($"  Receive         : {Poly2Bin(received, fq.Q)}");
            Console.WriteLine($"  Error Expected  : {Poly2Bin(err0, fq.Q)}");
            Console.WriteLine($"  Error Direct    : {Poly2Bin(err1, fq.Q)}");
            Console.WriteLine($"  Error Indirect  : {Poly2Bin(err2, fq.Q)}");
            Console.WriteLine($"  Nb Errors       : {err0.Coefs.Sum(a0 => a0.K)}");

            Console.WriteLine();
            if (!err2.Equals(err0) || !err1.Equals(err0))
                throw new();
        }
    }
    /// <summary>
    /// Runs a series of examples related to cyclotomic classes and BCH codes.
    /// </summary>
    public static void CyclotomicClassesExamples()
    {
        /*
         * Runs a series of examples related to cyclotomic classes and BCH codes.
         * The method demonstrates various use cases, including GAP examples, Wikipedia examples,
         * examples from "Algebre Tome1", and other miscellaneous examples.
         */

        // Set the polynomial display mode to use the 'StarCaret' format.
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;

        // GAP examples
        // Demonstrates the use of CyclotomicClasses method with different parameters (q, n).
        CyclotomicClasses(2, 21).Println(si => $"[{si.Glue("; ")}]", $"CycloCls (q:2, n:21)");
        CyclotomicClasses(10, 21).Println(si => $"[{si.Glue("; ")}]", $"CycloCls (q:10, n:21)");
        Console.WriteLine();

        // Wikipedia examples
        // Demonstrates the use of GeneratorBCHmδ method with different parameters (m, delta).
        GeneratorBCHmδ(4, 4);
        GeneratorBCHmδ(4, 6);
        GeneratorBCHmδ(4, 8);

        // Algebre Tome1, page 391-394
        // Demonstrates the use of GeneratorBCHmδ method with different parameters (m, delta).
        GeneratorBCHmδ(5, 7);
        GeneratorBCHmδ(7, 19);

        // Other examples
        // Demonstrates the use of GeneratorBCHmδ method with different parameters (m, delta).
        GeneratorBCHmδ(6, 5);
        GeneratorBCHmδ(6, 8);
    }

    /// <summary>
    /// Runs a series of detailed test examples for BCH codes with different parameters.
    /// </summary>
    public static void BCH_Codes_Examples()
    {
        /*
         * Runs a series of detailed test examples for BCH codes with different parameters.
         * The method includes multiple test cases using the CheckDetailedBCHmδ method
         * to simulate encoding, adding noise, and decoding processes.
         * Each test case specifies values for 'm', 'delta', and 'nbTest' parameters.
         */

        // Set the polynomial display mode to use the 'StarCaret' format.
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;

        // Test case 1: m = 5, delta = 7, nbTest = 10
        // This test provides a detailed example of encoding, adding noise,
        // and decoding for a BCH code with parameters (m: 5, delta: 7).
        // The word length is 17 bits, encoder length is 32 bits, and maximum errors allowed is 2 bits.
        CheckDetailedBCHmδ(m: 5, delta: 7, nbTest: 10);

        // Test case 2: m = 5, delta = 9, nbTest = 10
        // This test provides a detailed example of encoding, adding noise,
        // and decoding for a BCH code with parameters (m: 5, delta: 9).
        // The word length is 12 bits, encoder length is 32 bits, and maximum errors allowed is 3 bits.
        CheckDetailedBCHmδ(m: 5, delta: 9, nbTest: 10);

        // Test case 3: m = 7, delta = 19, nbTest = 10
        // This test provides a detailed example of encoding, adding noise,
        // and decoding for a BCH code with parameters (m: 7, delta: 19).
        // The word length is 72 bits, encoder length is 128 bits, and maximum errors allowed is 8 bits.
        CheckDetailedBCHmδ(m: 7, delta: 19, nbTest: 10);

        // Test case 4: m = 7, delta = 25, nbTest = 10
        // This test provides a detailed example of encoding, adding noise,
        // and decoding for a BCH code with parameters (m: 7, delta: 25).
        // The word length is 51 bits, encoder length is 128 bits, and maximum errors allowed is 11 bits.
        CheckDetailedBCHmδ(m: 7, delta: 25, nbTest: 10);
    }
    /*
        Syndrom : (a^5 + a^4 + a^2 + 1)*X^24 + (a^6 + a^5 + a^4 + a^3 + a^2)*X^23 + (a^5 + a^4 + a^3 + a^2 + 1)*X^22 + (a^3 + 1)*X^21 + (a^6 + a^4 + a^3)*X^20 + (a^6 + a^3 + a + 1)*X^19 + (a^5 + a^4 + a^3 + a^2 + 1)*X^18 + (a^3 + a^2 + a + 1)*X^17 + (a^5 + a^2 + a)*X^16 + (a^6 + a^5 + a^4 + a^2 + a + 1)*X^15 + (a^5 + a^4 + a^3 + a^2 + a)*X^14 + (a^6 + a + 1)*X^13 + (a^6 + a^3 + a^2 + a + 1)*X^12 + (a^6 + a^5 + a^3 + a + 1)*X^11 + (a^5 + a^3)*X^10 + (a^6 + a^5 + a^3 + a + 1)*X^9 + (a^6 + a^4 + a^3)*X^8 + (a^6 + a^5 + a^4 + a^3)*X^7 + (a^5 + a^4 + a^3 + a^2 + 1)*X^6 + (a^6 + a^5 + a^3 + a^2)*X^5 + (a^5 + a^3)*X^4 + (a^6 + a^5 + a^3 + a + 1)*X^3 + (a^6 + a^5 + a^3 + a^2)*X^2 + (a^6 + a^5 + a^2 + a)*X
          Direct          : u = X^9 + (a^5 + a^4)*X^8 + (a^4 + a)*X^7 + (a^6 + a^5 + a^4 + a^3 + a)*X^6 + (a^4 + a^2)*X^5 + (a^4 + a)*X^4 + (a^6 + a^5 + a^4 + a^2 + 1)*X^3 + (a^6 + a^5 + a^3 + a^2)*X^2 + a^3*X + a^6 + a^4 + a + 1 and v = X^9 + (a^4 + a)*X^7 + (a^4 + a^2)*X^5 + (a^6 + a^5 + a^4 + a^2 + 1)*X^3 + a^3*X
          BerlekampMassey : u = X^9 + (a^5 + a^4)*X^8 + (a^4 + a)*X^7 + (a^6 + a^5 + a^4 + a^3 + a)*X^6 + (a^4 + a^2)*X^5 + (a^4 + a)*X^4 + (a^6 + a^5 + a^4 + a^2 + 1)*X^3 + (a^6 + a^5 + a^3 + a^2)*X^2 + a^3*X + a^6 + a^4 + a + 1 and v = X^9 + (a^4 + a)*X^7 + (a^4 + a^2)*X^5 + (a^6 + a^5 + a^4 + a^2 + 1)*X^3 + a^3*X
          Verif Direct S  : True
          Verif B-M    S  : True
          diff u : 0
          diff v : 0
        Send Word8        : 11101000000110011110000011010110100100000010010111001000000001100010010000111011111110000000111100101111000000100011000011000010
          Receive         : 11101000000110011110000001110110100100000010010111011000000001100010010000011011111100010000101100101111000000110011001011000010
          Error Expected  : 00000000000000000000000010100000000000000000000000010000000000000000000000100000000010010000010000000000000000010000001000000000
          Error Direct    : 00000000000000000000000010100000000000000000000000010000000000000000000000100000000010010000010000000000000000010000001000000000
          Error Indirect  : 00000000000000000000000010100000000000000000000000010000000000000000000000100000000010010000010000000000000000010000001000000000
          Nb Errors       : 9

     */

    /// <summary>
    /// Runs a batch of tests for BCH codes with different parameters.
    /// </summary>
    public static void BCH_Codes_Batch()
    {
        /*
         * Runs a batch of tests for BCH codes with different parameters.
         * The method includes multiple test cases using the CheckBCHmδ method to simulate encoding, adding noise, and decoding processes.
         * Each test case specifies values for 'm', 'delta', 'nbTest', and 'maxErrors' parameters.
         */

        // Set the polynomial display mode to use the 'StarCaret' format.
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;

        // Test case 1: m = 5, delta = 7, nbTest = 100, maxErrors = true
        // This test simulates encoding, adding noise, and decoding for a BCH code with parameters (m: 5, delta: 7).
        // The word length is 17 bits, encoder length is 32 bits, and maximum errors allowed is 2 bits.
        CheckBCHmδ(m: 5, delta: 7, nbTest: 100, maxErrors: true);

        // Test case 2: m = 7, delta = 21, nbTest = 100, maxErrors = true
        // This test simulates encoding, adding noise, and decoding for a BCH code with parameters (m: 7, delta: 21).
        // The word length is 65 bits, encoder length is 128 bits, and maximum errors allowed is 9 bits.
        CheckBCHmδ(m: 7, delta: 21, nbTest: 100, maxErrors: true);
    }
    /*
    Word72     : 0011011110011110100100001110010110100110011111100110000110110100
      Send     : 00100100001000001100001010110101100110001011110011100001100011000010111010000010000000100101001111111110101101100101011001011000
      Received : 00100110001000101100101010110101100110001011111011110001100010001010111010000010000000100101011111111110101101100101011001111000
      Decoded  : 0011011110011110100100001110010110100110011111100110000110110100
      NbErrors : 9
      Verif    : True
     */
}