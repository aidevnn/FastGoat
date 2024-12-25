using System.Numerics;
using static FastGoat.Commons.IntExt;

namespace FastGoat.Commons;

/// <summary>
/// Provides extension methods for generating random distributions.
/// </summary>
public static class DistributionExt
{
    /// <summary>
    /// Generates a random integer between the specified minimum and maximum values, inclusive.
    /// </summary>
    /// <param name="min">The minimum value.</param>
    /// <param name="max">The maximum value.</param>
    /// <returns>A random integer between min and max.</returns>
    public static long Dice(long min, long max) => (long)double.Round((max - min) * Rng.NextDouble() + min);

    /// <summary>
    /// Generates a random integer between the specified minimum and maximum values, inclusive.
    /// </summary>
    /// <param name="min">The minimum value.</param>
    /// <param name="max">The maximum value.</param>
    /// <returns>A random integer between min and max.</returns>
    public static BigInteger Dice(BigInteger min, BigInteger max)
    {
        var v = (int)BigInteger.Log2(max - min) + 1;
        var rnd = v.SeqLazy().Select(_ => Rng.Next(2)).Aggregate(BigInteger.Zero, (acc, i) => 2 * acc + i);
        return BigInteger.Remainder(rnd, max - min + 1) + min;
    }

    /// <summary>
    /// Generates a sequence of random integers between the specified minimum and maximum values.
    /// </summary>
    /// <param name="size">The number of random integers to generate.</param>
    /// <param name="min">The minimum value.</param>
    /// <param name="max">The maximum value.</param>
    /// <returns>An IEnumerable of random integers.</returns>
    public static IEnumerable<long> DiceSample(int size, long min, long max) => size.SeqLazy().Select(_ => Dice(min, max));

    /// <summary>
    /// Generates a sequence of random integers between the specified minimum and maximum values.
    /// </summary>
    /// <param name="size">The number of random integers to generate.</param>
    /// <param name="min">The minimum value.</param>
    /// <param name="max">The maximum value.</param>
    /// <returns>An IEnumerable of random integers.</returns>
    public static IEnumerable<BigInteger> DiceSampleBigInt(int size, BigInteger min, BigInteger max) 
        => size.SeqLazy().Select(_ => Dice(min, max));

    /// <summary>
    /// Selects a random element from the specified array.
    /// </summary>
    /// <typeparam name="T">The type of the elements in the array.</typeparam>
    /// <param name="ts">The array to select from.</param>
    /// <returns>A random element from the array.</returns>
    public static T Dice<T>(T[] ts) => ts[Rng.Next(ts.Length)];
    
    /// <summary>
    /// Generates a sequence of random elements from the specified array.
    /// </summary>
    /// <typeparam name="T">The type of the elements in the array.</typeparam>
    /// <param name="size">The number of random elements to generate.</param>
    /// <param name="ts">The array to select from.</param>
    /// <returns>An IEnumerable of random elements.</returns>
    public static IEnumerable<T> DiceSample<T>(int size, T[] ts) => size.SeqLazy().Select(_ => Dice(ts));
    
    /// <summary>
    /// Performs a Bernoulli trial with the specified probability.
    /// </summary>
    /// <param name="p">The probability of success.</param>
    /// <returns>1 if the trial is successful; otherwise, 0.</returns>
    public static int Bernoulli(double p) => Rng.NextDouble() < p ? 1 : 0;
    
    /// <summary>
    /// Generates a sequence of Bernoulli trials with the specified probability.
    /// </summary>
    /// <param name="size">The number of trials to generate.</param>
    /// <param name="p">The probability of success.</param>
    /// <returns>An IEnumerable of trial results.</returns>
    public static IEnumerable<int> BernouilliSample(int size, double p) => size.SeqLazy().Select(_ => Bernoulli(p));
    
    /// <summary>
    /// Generates a random integer based on the Poisson distribution with the specified mean.
    /// </summary>
    /// <param name="l">The mean of the Poisson distribution.</param>
    /// <returns>A random integer based on the Poisson distribution.</returns>
    public static int Poisson(double l)
    {
        var L = Math.Exp(-l);
        int k = 0;
        double p = 1.0;
        while (p > L)
        {
            ++k;
            p *= Rng.NextDouble();
        }

        return k - 1;
    }
    
    /// <summary>
    /// Generates a sequence of random integers based on the Poisson distribution with the specified mean.
    /// </summary>
    /// <param name="size">The number of random integers to generate.</param>
    /// <param name="p">The mean of the Poisson distribution.</param>
    /// <returns>An IEnumerable of random integers.</returns>
    public static IEnumerable<int> PoissonSample(int size, double p) => size.SeqLazy().Select(_ => Poisson(p));

    /// <summary>
    /// Generates a random integer based on the binomial distribution with the specified parameters.
    /// </summary>
    /// <param name="n">The number of trials.</param>
    /// <param name="p">The probability of success.</param>
    /// <returns>A random integer based on the binomial distribution.</returns>
    public static int Binomial(int n, double p) => n.SeqLazy().Select(_ => Bernoulli(p)).Sum();
    
    /// <summary>
    /// Generates a sequence of random integers based on the binomial distribution with the specified parameters.
    /// </summary>
    /// <param name="size">The number of random integers to generate.</param>
    /// <param name="n">The number of trials.</param>
    /// <param name="p">The probability of success.</param>
    /// <returns>An IEnumerable of random integers.</returns>
    public static IEnumerable<int> BinomialSample(int size, int n, double p) 
        => size.SeqLazy().Select(_ => Binomial(n, p));

    /// <summary>
    /// Generates a centered binomial random integer with the specified parameters.
    /// </summary>
    /// <param name="n">The number of trials.</param>
    /// <param name="p">The probability of success.</param>
    /// <returns>A centered binomial random integer.</returns>
    public static int CenteredBinomial(int n, double p) => n.SeqLazy().Select(_ => Bernoulli(p)).Sum() - n / 2;
    
    /// <summary>
    /// Generates a sequence of centered binomial random integers with the specified parameters.
    /// </summary>
    /// <param name="size">The number of random integers to generate.</param>
    /// <param name="n">The number of trials.</param>
    /// <param name="p">The probability of success.</param>
    /// <returns>An IEnumerable of centered binomial random integers.</returns>
    public static IEnumerable<int> CenteredBinomialSample(int size, int n, double p) 
        => size.SeqLazy().Select(_ => CenteredBinomial(n, p));
    
    /// <summary>
    /// Generates a binomial random integer with equal probability of success and failure.
    /// </summary>
    /// <param name="n">The number of trials.</param>
    /// <returns>A binomial random integer.</returns>
    public static int BinomialEquiProb(int n) => Binomial(n, 0.5);
    
    /// <summary>
    /// Generates a sequence of binomial random integers with equal probability of success and failure.
    /// </summary>
    /// <param name="size">The number of random integers to generate.</param>
    /// <param name="n">The number of trials.</param>
    /// <returns>An IEnumerable of binomial random integers.</returns>
    public static IEnumerable<int> BinomialEquiProbSample(int size, int n) 
        => size.SeqLazy().Select(_ => BinomialEquiProb(n));
    
    /// <summary>
    /// Generates a centered binomial random integer with equal probability of success and failure.
    /// </summary>
    /// <param name="n">The number of trials.</param>
    /// <returns>A centered binomial random integer.</returns>
    public static int CenteredBinomialEquiProb(int n) => CenteredBinomial(n, 0.5);
    
    /// <summary>
    /// Generates a sequence of centered binomial random integers with equal probability of success and failure.
    /// </summary>
    /// <param name="size">The number of random integers to generate.</param>
    /// <param name="n">The number of trials.</param>
    /// <returns>An IEnumerable of centered binomial random integers.</returns>
    public static IEnumerable<int> CenteredBinomialEquiProbSample(int size, int n) 
        => size.SeqLazy().Select(_ => CenteredBinomialEquiProb(n));

    /// <summary>
    /// Generates a pair of random numbers based on the normal distribution using the Box-Muller transform.
    /// </summary>
    /// <param name="mu">The mean of the distribution.</param>
    /// <param name="sigma">The standard deviation of the distribution.</param>
    /// <returns>A tuple containing two random numbers based on the normal distribution.</returns>
    public static (double n1, double n2) NormalBoxMuller(double mu, double sigma)
    {
        var R = Double.Sqrt(-2 * Double.Log(Rng.NextDouble()));
        var (sin, cos) = Double.SinCosPi(2 * Rng.NextDouble());
        return (R * sigma * sin + mu, R * sigma * cos + mu);
    }
    
    /// <summary>
    /// Generates a pair of random numbers based on the standard normal distribution using the Box-Muller transform.
    /// </summary>
    /// <returns>A tuple containing two random numbers based on the standard normal distribution.</returns>
    public static (double n1, double n2) StdNormalBoxMuller() => NormalBoxMuller(0.0, 1.0);
    
    /// <summary>
    /// Generates a sequence of random numbers based on the normal distribution using the Box-Muller transform.
    /// </summary>
    /// <param name="size">The number of random numbers to generate.</param>
    /// <param name="mu">The mean of the distribution.</param>
    /// <param name="sigma">The standard deviation of the distribution.</param>
    /// <returns>An IEnumerable of random numbers based on the normal distribution.</returns>
    public static IEnumerable<double> NormalSample(int size, double mu, double sigma)
    {
        var half = size / 2 - 1;
        var isPair = size % 2 == 0;
        for (int i = 0; i <= half; i++)
        {
            var (n1, n2) = NormalBoxMuller(mu, sigma);
            yield return n1;
            if (isPair || i < half)
                yield return n2;
        }
    }

    /// <summary>
    /// Generates a sequence of random numbers based on the standard normal distribution using the Box-Muller transform.
    /// </summary>
    /// <param name="size">The number of random numbers to generate.</param>
    /// <returns>An IEnumerable of random numbers based on the standard normal distribution.</returns>
    public static IEnumerable<double> StdNormalSample(int size) => NormalSample(size, 0.0, 1.0);

    public static int DiscreteGaussian(double mu, double sigma, double tau)
    {
        var (n0, n1) = NormalBoxMuller(mu, sigma);
        if (double.Abs(n0 - mu) < tau * sigma)
            return (int)double.Round(n0);
        if (double.Abs(n1 - mu) < tau * sigma)
            return (int)double.Round(n1);

        return 0;
    }

    public static IEnumerable<int> DiscreteGaussianSample(int size, double mu, double sigma, double tau)
        => size.SeqLazy().Select(_ => DiscreteGaussian(mu, sigma, tau));

    public static IEnumerable<int> DiscreteGaussianSample(int size, double mu, double sigma) =>
        DiscreteGaussianSample(size, mu, sigma, size / sigma);
}