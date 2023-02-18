using System.Diagnostics;

namespace FastGoat.Commons;

public static class GlobalStopWatch
{
    static GlobalStopWatch()
    {
        sw = Stopwatch.StartNew();
    }

    static Stopwatch sw { get; }
    public static void Restart() => sw.Restart();
    public static void Stop() => sw.Stop();
    public static void Show(string label) => Console.WriteLine($"# {label} Time:{sw.ElapsedMilliseconds} ms");

    public static long Time(string label, Action action, bool lbl = true)
    {
        if (lbl)
            Console.WriteLine($"# {label} Start");
        
        sw.Restart();
        action();
        sw.Stop();
        Console.WriteLine($"# {label} Time:{sw.ElapsedMilliseconds} ms");
        return sw.ElapsedMilliseconds;
    }

    public static void Bench(int nb, string label, Action action)
    {
        if (nb < 1)
            throw new();
        
        var list = new List<long>();
        for (int i = 0; i < nb; i++)
            list.Add(Time(label, action, false));

        var avg = list.Average();
        var dev = Double.Sqrt(list.Select(t => Double.Pow(t - avg, 2)).Average());
        Console.WriteLine($"# {label} Avg Time:{(long)avg} ms Dev:{dev:F}");
        Console.WriteLine();
    }

    private static int ct = 0;
    public static void InfiniteLoopBreakerReset() => ct = 0;

    public static void InfiniteLoopBreaker(int n, string msg = "")
    {
        ++ct;
        if (ct > n)
            throw new ArgumentException($"################## Infinite Loop breaker ### {msg}");
    }
}